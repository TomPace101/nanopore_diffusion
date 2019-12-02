"""Solve the multi-species unhomogenized Fickian diffusion problem,
and extract the data needed for post-processing efforts.
A single species is assumed, as there is no interaction or potential.
Isotropy is assumed, but the diffusion constant may vary spatially."""

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
from . import common_methods

class FLSimulator(simrequest.SimulationRequest):
  """Simulator for for Unhomogenized Fickian Diffusion"""

  #Common methods
  calcflux = common_methods.calcflux
  fluxintegral = common_methods.fluxintegral

  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.meshinfo.mesh,'CG',conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.meshinfo.mesh, "CG", conditions.elementorder)

    #Trial Function
    self.c = fem.TrialFunction(self.V)

    #Test Function
    v=fem.TestFunction(self.V)

    #Solution function
    self.soln=fem.Function(self.V,name="soln")

    #Measure and normal for external boundaries
    self.ds = fem.Measure("ds",domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n=fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Load diffusion coefficient Function
    self.Dlocal=fem.Function(self.V,name="Dlocal")
    self.process_load_commands()

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]

    #Neumann boundary conditions
    self.nbcs = {}
    neumann=getattr(conditions,'neumann',{})
    for psurf,value in neumann.items():
      if value is not None:
        if type(value)==int or type(value)==float:
          if value != 0: #Neumann conditions of zero can be omitted from the weak form, to the same effect
            self.nbcs[psurf]=fem.Constant(value)
        elif type(value)==list:
          exprstr, exprargs = value
          self.nbcs[psurf]=fem.Expression(exprstr,element=ele,**exprargs)
    
    #Weak Form
    allterms=equationbuilder.EquationTermDict()
    #Body term
    termname='body'
    allterms.add(termname,-self.Dlocal*fem.dot(fem.grad(self.c),fem.grad(v))*self.dx,bilinear=True)
    #If no boundary terms will be added, go ahead and apply zero
    if len(self.nbcs)==0:
      termname='boundary_zero'
      allterms.add(termname,fem.Constant(0)*v*self.dx,bilinear=False)
    #Boundary terms for Neumann conditions
    for psurf,expr in self.nbcs.items():
      termname='boundary_neumann_%d'%(psurf)
      bterm = self.Dlocal*expr*v*self.ds(psurf)
      allterms.add(termname,bterm,bilinear=False)

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

  def effective_D(self,outattr,fluxattr,areaattr,startloc,endloc,solnattr='soln',idx=None):
    """Calculate effective diffusion constant

    Arguments:

      - outattr = attribute path for storage of result
      - fluxattr = attribute path to previously calculated total

          This requires a previous call to fluxintegral.

      - areaattr = attribute path to previously calculated area in results dictionary

          This requires a previous call to facet_area.

      - startloc = argument to get_pointcoords for start of line
      - endloc = argument to get_pointcoords for end of line
      - solnattr = optional, attribute path to the concentration solution
      - idx = index of the solution field to write out, None (default) if not a sequence

    No return value.
    No output files."""
    #Get the flux and area values
    totflux=self.get_nested(fluxattr)
    area=self.get_nested(areaattr)
    #Get the object with the solution data
    vals=self.get_nested(solnattr)
    if idx is not None:
      vals = vals[idx]
    #Get the points for data extraction
    assert len(startloc)==len(endloc), "Start and end locations have different dimensionality"
    startcoords=self.get_pointcoords(startloc)
    endcoords=self.get_pointcoords(endloc)
    #Calculate distance between the two points
    deltas=[p[1]-p[0] for p in zip(startcoords,endcoords)]
    delta_s=np.sqrt(sum([d**2 for d in deltas]))
    #Calculate the change in concentration between the two points
    delta_c=vals(*endcoords)-vals(*startcoords)
    #Calculate diffusion constant
    Deff=float(totflux/area*delta_s/delta_c)
    #Store result
    self.set_nested(outattr,Deff)
    return

#Register for loading from yaml
yaml_manager.register_classes([FLSimulator])




