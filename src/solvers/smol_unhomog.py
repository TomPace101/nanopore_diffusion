#Solve the unhomogenized Smoluchowski diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import math
import os
import os.path as osp

#Site packages
import fenics as fem

#Local
import solver_general

class LPBConditions(solver_general.GenericConditions):
  """Condition defnitions for use with LPBSolver
  Attributes:
    dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    debye_length = Debye length"""
  __slots__=['dirichlet','debye_length']

class LPBSolver(solver_general.GenericSolver):
  """Solver for linearized Poisson-Boltzmann equation.
  Additional attributes not inherited from GenericSolver:
    conditions = instance of LPBConditions
    lambda_D = Debye length
    V = FEniCS FunctionSpace on the mesh
    bcs = list of FEniCS DirichletBC instances
    ds = FEniCS Measure for facet boundary conditions
    phi = FEniCS TrialFunction on V
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams,other):
    """Initialize the model.
    Arguments:
      modelparams = solver_run.ModelParameters instance
      meshparams = buildgeom.MeshParameters instance
      other = solver to get mesh from"""

    #Load parameters, init output, mesh setup
    super(LPBSolver, self).__init__(modelparams,meshparams)

    #Load mesh and meshfunctions
    self.mesh=other.mesh
    self.facets=other.facets
    self.cells=other.cells
    self.V = other.V
    self.ds = other.ds

    #Get conditions
    self.conditions=LPBConditions(**modelparams.conditions)

    #Properties of problem domain
    self.lambda_D = self.conditions.debye_length

    #Function space for scalars and vectors
    ##self.V = fem.FunctionSpace(self.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.facets,psurf) for psurf,val in self.conditions.dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    ##self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.facets) ##TODO: specify which facets are Neumann?
    if hasattr(self.conditions,'neumann'):
      raise NotImplementedError

    #Define variational problem
    self.phi=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    self.a=((1/self.lambda_D**2)*self.phi*self.v + fem.dot(fem.grad(self.phi),fem.grad(self.v)))*fem.dx
    self.L=fem.Constant(0)*self.v*self.ds

  def solve(self):
    "Do the step of solving this equation"
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bcs)
    return

#Lookup of electric potential solvers by name
potentialsolverclasses={'linear_pb':LPBSolver}

class SUConditions(solver_general.GenericConditions):
  """Condition defnitions for use with SUSolver
  Attributes:
    dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    D_bulk = bulk diffusion constant
    q = electric charge of ion
    beta = 1/kBT for the temperature under consideration, in units compatible with q times the potential
    potential = dictionary defining solver_run.ModelParameters for electric potential
    trans_dirichlet = Dirichlet boundary conditions after Slotboom transformation
  Note also that the attribute bclist (inherited), contains Dirichlet conditions on c, rather than cbar.
    That is, the code will do the Slotboom transformation on the Dirichlet boundary conditions."""
  __slots__=['dirichlet','D_bulk','q','beta','potential','trans_dirichlet']
  def transform_bcs(self,pdict,beta_q):
    """Apply Slotboom transformation to Dirichlet boundary conditions.
    This function requires that the facets with Dirichlet conditions for the concentration
      must also have Dirichlet conditions for the potential. (The reverse is not required.)
    Arguments:
      pdict = dictionary of Dirichlet boundary conditions for the potential
      beta_q = product of beta and q
    No return value.
    trans_dirichlet attribute is updated"""
    transvals={}
    for psurf,cval in self.dirichlet.items():
      transvals[psurf] = cval*math.exp(beta_q*pdict[psurf])
    self.trans_dirichlet=transvals
    return

class SUSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized Smoluchowski Diffusion
  Additional attributes not inherited from GenericSolver:
    conditions = instance of SUConditions
    beta_q = product of beta and q (for convenience)
    V = FEniCS FunctionSpace on the mesh
    V_vec = FEniCS VectorFunctionSpace on the mesh
    bcs = FEniCS BCParameters
    ds = FEniCS Measure for facet boundary conditions
    potsolv = instance of solver for the electric potential
    Dbar = FEniCS Function
    cbar = FEniCS TrialFunction on V
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams):
    """Initialize the model.
    Arguments:
      modelparams = solver_run.ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""

    #Load parameters, init output, mesh setup
    super(SUSolver, self).__init__(modelparams,meshparams)
    self.loadmesh()

    #Get conditions
    self.conditions=SUConditions(**modelparams.conditions)
    self.beta_q = self.conditions.beta * self.conditions.q

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.mesh, "CG", self.conditions.elementorder)

    #Measure for external boundaries
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.facets)

    #Set up electric potential field
    potentialparams_dict=self.conditions.potential
    for key in ['modelname','meshname','meshparamsfile','basename']:
      potentialparams_dict[key]=getattr(modelparams,key)
    potentialparams=solver_general.ModelParametersBase(**potentialparams_dict)
    self.potsolv=potentialsolverclasses[potentialparams.equation](potentialparams,meshparams,self)
    self.potsolv.diskwrite=False
    self.potsolv.solve()
    self.potsolv.create_output()
    self.info['potential']=self.potsolv.info
    self.outdata.plots=self.potsolv.outdata.plots

    #Dirichlet boundary conditions
    self.conditions.transform_bcs(self.potsolv.conditions.dirichlet,self.beta_q) #apply Slotboom transformation
    self.bcs=[fem.DirichletBC(self.V,val,self.facets,psurf) for psurf,val in self.conditions.trans_dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    if hasattr(self.conditions,'neumann'):
      raise NotImplementedError

    #Define the Dbar function
    self.Dbar=self.conditions.D_bulk*fem.exp(-self.beta_q*self.potsolv.soln)
    self.Dbar_proj=fem.project(self.Dbar,self.V)

    #Define variational problem
    self.d3x = fem.Measure('cell',domain=self.mesh)
    self.cbar=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    self.a=self.Dbar*fem.dot(fem.grad(self.cbar),fem.grad(self.v))*self.d3x
    self.L=fem.Constant(0)*self.v*self.ds

  def solve(self):
    "Do the step of solving this equation"
    #Solve for cbar
    self.sb_soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.sb_soln, self.bcs)
    #Transform back to concentration
    self.soln=fem.project(self.sb_soln*fem.exp(-self.beta_q*self.potsolv.soln),self.V)
    return

solverclasses={'smol_unhomog':SUSolver}
