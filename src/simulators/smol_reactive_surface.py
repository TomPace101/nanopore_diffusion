"""Solve the unhomogenized Smoluchowski diffusion problem
and extract the data needed for post-processing efforts"""

#Standard library
import argparse
import math
import os
import os.path as osp

#Site packages
import fenics as fem

#Local
import common
import simulator_general

class LPBConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with LPBSimulator

  Attributes:

    - dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    - debye_length = Debye length"""
  __slots__=['dirichlet','debye_length']

class LPBSimulator(simulator_general.GenericSimulator):
  """Simulator for linearized Poisson-Boltzmann equation.

  Additional attributes not inherited from GenericSimulator:

    - conditions = instance of LPBConditions
    - lambda_D = Debye length
    - V = FEniCS FunctionSpace on the mesh
    - bcs = list of FEniCS DirichletBC instances
    - ds = FEniCS Measure for facet boundary conditions
    - phi = FEniCS TrialFunction on V
    - v = FEniCS TestFunction on V
    - a = bilinear form in variational problem
    - L = linear form in variational problem"""
  def __init__(self,modelparams,other):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance
      - other = simulator to get mesh from"""

    #Load parameters, init output, mesh setup
    super(LPBSimulator, self).__init__(modelparams)

    #Load mesh and meshfunctions
    self.meshinfo=other.meshinfo
    self.V = other.V_scalar
    self.ds = other.ds

    #Get conditions
    self.conditions=LPBConditions(**modelparams.conditions)

    #Properties of problem domain
    self.lambda_D = self.conditions.debye_length

    #Function space for scalars and vectors
    ##self.V = fem.FunctionSpace(self.meshinfo.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in self.conditions.dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    ##self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets) ##TODO: specify which facets are Neumann?
    if hasattr(self.conditions,'neumann'):
      raise NotImplementedError

    #Define variational problem
    self.phi=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    self.a=((1/self.lambda_D**2)*self.phi*self.v + fem.dot(fem.grad(self.phi),fem.grad(self.v)))*fem.dx
    self.L=fem.Constant(0)*self.v*self.ds

  def run(self):
    "Run this simulation"
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bcs)
    return

#Lookup of electric potential simulators by name
potentialsimulatorclasses={'linear_pb':LPBSimulator}

class EquationTerm(simulator_general.EquationTermBase):
  __slots__=('bilinear','termid','species','bound_surf')


class Species(common.ParameterSet):
  """Information for a single chemical species.
  
  User-Provided attributes:
  
    - symbol: chemical symbol, as string
    - z: ionic charge, as number
    - D: diffusion constant, as number"""
  __slots__=('symbol','z','D')
  _required_attrs=__slots__

class SUConditions(simulator_general.GenericConditions):
  """Condition defnitions for use with SUSimulator

  User-Provided Attributes:

    - dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    - species = sequence of dictionaries, each defining a Species object
    - beta = 1/kBT for the temperature under consideration, in units compatible with q times the potential
    - potential = dictionary defining simulator_run.ModelParameters for electric potential
  
  Note also that the attribute bclist (inherited), contains Dirichlet conditions on c, rather than cbar.
    That is, the code will do the Slotboom transformation on the Dirichlet boundary conditions."""
  __slots__=['dirichlet','species','beta','potential','trans_dirichlet']

class SUSimulator(simulator_general.GenericSimulator):
  """Simulator for Unhomogenized Smoluchowski Diffusion

  Additional attributes not inherited from GenericSimulator:

    - conditions = instance of SUConditions
    - V = FEniCS FunctionSpace on the mesh
    - V_vec = FEniCS VectorFunctionSpace on the mesh
    - bcs = FEniCS BCParameters
    - ds = FEniCS Measure for facet boundary conditions
    - potsim = instance of simulator for the electric potential
    - Dbar = FEniCS Function
    - cbar = FEniCS TrialFunction on V
    - v = FEniCS TestFunction on V
    - a = bilinear form in variational problem
    - L = linear form in variational problem"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(SUSimulator, self).__init__(modelparams)

    #Get conditions
    self.conditions=SUConditions(**modelparams.conditions)

    #Species
    self.species=[]
    for s,d in enumerate(self.conditions.species):
      spec=Species(**d)
      self.species.append(spec)
    self.Nspecies=len(self.species)
    
    #Elements and Function space(s)
    ele = fem.FiniteElement('CG',self.meshinfo.mesh.ufl_cell(),self.conditions.elementorder)
    mele = fem.MixedElement([ele]*self.Nspecies)
    self.V = fem.FunctionSpace(self.meshinfo.mesh,mele)
    self.V_scalar=fem.FunctionSpace(self.meshinfo.mesh,'CG',self.conditions.elementorder)

    #Trial Functions
    self.u = fem.TrialFunction(self.V)
    cbarlist=fem.split(self.u)

    #Test Functions
    vlist=fem.TestFunctions(self.V)

    #Solution function(s)
    self.soln=fem.Function(self.V)

    #Measure and normal for external boundaries
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n=fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Set up electric potential field
    potentialparams_dict=self.conditions.potential
    for key in ['modelname','meshname','basename']:
      potentialparams_dict[key]=getattr(modelparams,key)
    potentialparams=simulator_general.ModelParametersBase(**potentialparams_dict)
    self.potsim=potentialsimulatorclasses[potentialparams.equation](potentialparams,self)
    self.potsim.diskwrite=False
    self.potsim.run()
    self.potsim.create_output()
    self.info['potential']=self.potsim.info
    self.outdata.plots=self.potsim.outdata.plots

    #Slootboom transformation for dirichlet condition
    def transform_value(cval,phival,beta_z):
      return cval*math.exp(beta_z*phival)

    #Dirichlet boundary conditions
    pot_d=self.potsim.conditions.dirichlet
    self.bcs=[]
    dirichlet=getattr(self.conditions,'dirichlet',{})
    for psurf,vals in dirichlet.items():
      for s,value in enumerate(vals):
        if value is not None:
          transval=transform_value(value,pot_d[psurf],self.conditions.beta*self.species[s].z)
          self.bcs.append(fem.DirichletBC(self.V.sub(s),fem.Constant(transval),self.meshinfo.facets,psurf))

    #Neumann boundary conditions
    #assumed to all be zero in this case
    if hasattr(self.conditions,'neumann'):
      raise NotImplementedError

    #Calculate Dbar for each species
    self.Dbar_dict={}
    self.Dbar_proj=[]
    for s,spec in enumerate(self.species):
      Dbar=spec.D*fem.exp(-self.conditions.beta*spec.z*self.potsim.soln)
      self.Dbar_dict[s]=Dbar
      self.Dbar_proj.append(fem.project(Dbar,self.V_scalar))

    #Weak Form
    allterms=simulator_general.EquationTermDict(EquationTerm)
    for s,cbar in enumerate(cbarlist):
      if self.species[s].D is not None:
        termname='body_%d'%s
        allterms.add(termname,self.Dbar_dict[s]*fem.dot(fem.grad(cbar),fem.grad(vlist[s]))*self.dx,bilinear=True)

        termname='boundary_%d'%s
        allterms.add(termname,fem.Constant(0)*vlist[s]*self.dx,bilinear=False)

    #Problem and Solver
    self.a=allterms.sumterms(bilinear=True)
    self.L=allterms.sumterms(,bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

  def run(self):
    "Run this simulation"
    #solve
    self.solver.solve()
    #transform back
    self.solnlist=fem.split(self.soln)
    self.clist=[]
    for s,cbar in enumerate(self.solnlist):
      c=fem.project(cbar*fem.exp(-self.conditions.beta*self.species[s].z*self.potsim.soln),self.V_scalar)
      self.clist.append(c)
    return

simulatorclasses={'smol_reactive_surface':SUSimulator}
