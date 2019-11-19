"""Solve the multi-species unhomogenized Smoluchowski diffusion problem,
with a reactive boundary condition,
and extract the data needed for post-processing efforts.

Species may vary within each subdomain."""

#Standard library
from argparse import Namespace
from collections import OrderedDict as odict
import math

#Site packages
import numpy as np
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder
from . import common_methods

_SUConditions_props_schema_yaml="""#SUConditions
species:
  type: array
  items: {type: object}
reactive:
  type: object
  items:
    type: array
    items:
      - {type: string}
      - {type: string}
beta: {type: number}
potential_dirichlet: {type: object}
"""
SUConditions_props_schema=yaml_manager.readstring(_SUConditions_props_schema_yaml)
SUConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    SUConditions_props_schema,['species'])

class SUSimulator(simrequest.SimulationRequest):
  """Simulator for for Unhomogenized Smoluchowski Diffusion"""
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,SUConditions_schema)

  #Common methods
  calcflux = common_methods.calcflux
  fluxintegral = common_methods.fluxintegral

  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Species
    self.species=[]
    self.species_dict=odict()
    self.species_indices=odict()
    for s,d in enumerate(conditions.species):
      spec=Namespace(**d)
      self.species.append(spec)
      self.species_dict[spec.symbol]=spec
      self.species_indices[spec.symbol]=s
    self.Nspecies=len(self.species)
    
    #Elements and Function space(s)
    ele = fem.FiniteElement('CG',self.meshinfo.mesh.ufl_cell(),conditions.elementorder)
    mele = fem.MixedElement([ele]*self.Nspecies)
    self.V = fem.FunctionSpace(self.meshinfo.mesh,mele)
    self.V_scalar=fem.FunctionSpace(self.meshinfo.mesh,'CG',conditions.elementorder)
    self.V_vec=fem.VectorFunctionSpace(self.meshinfo.mesh,'CG',conditions.elementorder)

    #Trial Functions
    self.u = fem.TrialFunction(self.V)
    cbarlist=fem.split(self.u)

    #Test Functions
    vlist=fem.TestFunctions(self.V)

    #Solution function(s)
    self.soln=fem.Function(self.V,name="soln")

    #Measure and normal for external boundaries
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n=fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Load electric potential as a Function
    self.potential=fem.Function(self.V_scalar,name="Phi")
    self.process_load_commands()

    #Slootboom transformation for dirichlet condition
    def transform_value(cval,phival,beta_z):
      return cval*math.exp(beta_z*phival)

    #Dirichlet boundary conditions
    ##TODO: this requires any dirichlet boundary condition for concentration to also have a dirichlet boundary condition for the potential
    #that way, there is a constant dirichlet boundary condition for cbar
    pot_d=conditions.potential_dirichlet
    self.bcs=[]
    dirichlet=getattr(conditions,'dirichlet',{})
    for psurf,vals in dirichlet.items():
      for s,value in enumerate(vals):
        if value is not None:
          transval=transform_value(value,pot_d[psurf],conditions.beta*self.species[s].z)
          fspace=self.V.sub(s) if self.Nspecies>1 else self.V
          setvalue=transval if self.Nspecies>1 else (transval,)
          self.bcs.append(fem.DirichletBC(fspace,fem.Constant(setvalue),self.meshinfo.facets,psurf))

    #Neumann boundary conditions
    self.nbcs = {}
    neumann=getattr(conditions,'neumann',{})
    for psurf,vals in neumann.items():
      for s,value in enumerate(vals):
        if value is not None:
          if type(value)==int or type(value)==float:
            if value != 0: #Neumann conditions of zero can be omitted from the weak form, to the same effect
              self.nbcs[(psurf,s)]=fem.Constant(value)
          elif type(value)==list:
            exprstr, exprargs = value
            self.nbcs[(psurf,s)]=fem.Expression(exprstr,element=ele,**exprargs)
    
    #Reactive boundary terms
    #just convert symbols to species indices
    self.rbcs = {}
    reactive=getattr(conditions,'reactive',{})
    for psurf,pair in reactive.items():
      spair=[self.species_indices[symb] for symb in pair]
      self.rbcs[psurf]=spair

    #Calculate Dbar for each species
    self.Dbar_dict={}
    self.Dbar_proj=[]
    for s,spec in enumerate(self.species):
      Dbar=spec.D*fem.exp(-conditions.beta*spec.z*self.potential)
      self.Dbar_dict[s]=Dbar
      thisproj=fem.project(Dbar,self.V_scalar,solver_type="cg",preconditioner_type="amg") #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)
      thisproj.rename('Dbar','transformed diffusion coefficient')
      self.Dbar_proj.append(thisproj)

    #Weak Form
    allterms=equationbuilder.EquationTermDict()
    #Body terms
    for s,cbar in enumerate(cbarlist):
      if self.species[s].D is not None:
        termname='body_%d'%s
        allterms.add(termname,self.Dbar_dict[s]*fem.dot(fem.grad(cbar),fem.grad(vlist[s]))*self.dx,bilinear=True)
        #Not knowing if any boundary term will be added here, go ahead and apply zero
        termname='boundary_zero_%d'%s
        allterms.add(termname,fem.Constant(0)*vlist[s]*self.dx,bilinear=False)
    #Boundary terms for Neumann conditions
    for tup,expr in self.nbcs.items():
      psurf,s = tup
      termname='boundary_neumann_%d_%d'%(s,psurf)
      bterm = expr*vlist[s]*self.ds(psurf)
      allterms.add(termname,bterm,bilinear=False)
    #Boundary terms for reactive boundaries
    for psurf,pair in self.rbcs.items():
      r,p=pair
      termname='reactive_%d'%psurf
      bterm = self.Dbar_dict[r]*fem.dot(self.n,fem.grad(cbarlist[r]))*(vlist[r]-vlist[p])*self.ds(psurf)
      allterms.add(termname,-bterm,bilinear=True)

    #Problem and Solver
    self.a=allterms.sumterms(bilinear=True)
    self.L=allterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)
    
    #Get solver parameters
    self.set_solver_parameters()

    #solve
    if not getattr(self,'skipsolve',False):
      self.solver.solve()

    #transform back
    self.solnlist=fem.split(self.soln)
    self.clist=[]
    self.cbarlist=[]
    for s,cbar in enumerate(self.solnlist):
      symb=self.species[s].symbol
      expr=cbar*fem.exp(-conditions.beta*self.species[s].z*self.potential)
      c=fem.project(expr,self.V_scalar,solver_type="cg",preconditioner_type="amg")
      c.rename('c_'+symb,'concentration of species '+symb)
      self.clist.append(c)
      cbar_single=fem.project(cbar,self.V_scalar,solver_type="cg",preconditioner_type="amg")
      cbar_single.rename('cbar_'+symb,'transformed concentration of species '+symb)
      self.cbarlist.append(cbar_single)



#Register for loading from yaml
yaml_manager.register_classes([SUSimulator])

