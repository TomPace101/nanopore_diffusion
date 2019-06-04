"""Solve the linearized Poisson-Boltzmann equation"""

#Standard library
from argparse import Namespace

#Site packages
import numpy as np
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder

_LPBConditions_props_schema_yaml="""#LPBConditions
kappa: {type: number}
"""
LPBConditions_props_schema=yaml_manager.readstring(_LPBConditions_props_schema_yaml)
LPBConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    LPBConditions_props_schema,['kappa'])

class LPBSimulator(simrequest.SimulationRequest):
  """Simulator for linearized Poisson-Boltzmann equation
  
  User-defined attributes:
  
    - """
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,LPBConditions_schema)    
  
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Properties of problem domain
    self.kappa = self.conditions.kappa
    
    #Function Spaces and Functions
    #Function spaces
    self.V=fem.FunctionSpace(self.meshinfo.mesh,'CG',conditions.elementorder)
    #Trial and Test Functions
    self.phi=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V)

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    if hasattr(conditions,'neumann'):
      raise NotImplementedError

    #Measures and facet normals
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)

    #Equation term dictionary
    self.eqnterms=equationbuilder.EquationTermDict()

    #Weak form
    #kappa term
    ufl_term1=(self.kappa**2)*self.phi*self.v*fem.dx
    if self.kappa is not None:
      self.eqnterms.add('bilinear_1',ufl_term1,bilinear=True)
    #Second term
    ufl_term2=fem.dot(fem.grad(self.phi),fem.grad(self.v))*fem.dx
    self.eqnterms.add('bilinear_2',ufl_term2,bilinear=True)
    #Linear forms
    self.eqnterms.add('linear',fem.Constant(0)*self.v*self.ds,bilinear=False)

    #Problem and Solver
    self.a=self.eqnterms.sumterms(bilinear=True)
    self.L=self.eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

    #Get solver parameters from the conditions
    for k,v in getattr(conditions,'solver_parameters',{}).items():
      self.solver.parameters[k]=v

    #Solve
    self.solver.solve()

#Register for loading from yaml
yaml_manager.register_classes([LPBSimulator])
