"""Solve the homogenized Fickian diffusion problem (corrector problem)
and extract the data needed for post-processing efforts"""

#Standard library
from argparse import Namespace

#Site packages
import numpy as np
import fenics as fem
import ufl

#This package
from ..requesthandler import yaml_manager
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder

BOUNDTOL=1e-6

class PeriodicBoundary2D(fem.SubDomain):
  """SubDomain subclass for 2D Periodic boundary condition"""
  def __init__(self,xlims,ylims):
    """Arguments:
    
      - xlims = pair of x-values: (xmin,xmax)
      - ylims = pair of y-values: (ymin,ymax)"""
    super(PeriodicBoundary2D, self).__init__()
    self.left,  self.right = xlims
    self.bottom,self.top   = ylims
    self.xspan = self.right-self.left
    self.yspan = self.top-self.bottom
  # Left boundary is "target domain" G
  def inside(self, x, on_boundary):
    # return True if on left or bottom boundary AND NOT on one of the two corners 
    #  (self.left, self.top) and (self.right, self.bottom)
    return bool((fem.near(x[0], self.left) or fem.near(x[1], self.bottom)) and 
                (not ((fem.near(x[0], self.left) and fem.near(x[1], self.top)) or 
                      (fem.near(x[0], self.right) and fem.near(x[1], self.bottom)))) and on_boundary)
  def map(self, x, y):
    if fem.near(x[0], self.right) and fem.near(x[1], self.top):
      y[0] = x[0] - self.xspan
      y[1] = x[1] - self.yspan
    elif fem.near(x[0], self.right):
      y[0] = x[0] - self.xspan
      y[1] = x[1]
    else:   # fem.near(x[1], self.top)
      y[0] = x[0]
      y[1] = x[1] - self.yspan

##TODO
#Get the class below working
# class PeriodicBoundary3D(fem.SubDomain):
#   """SubDomain subclass for 3D Periodic boundary condition"""
#   def __init__(self,xlims,ylims,zlims):
#     """Arguments:
# 
#       - xlims = pair of x-values: (xmin,xmax)
#       - ylims = pair of y-values: (ymin,ymax)
#       - zlims = pair of z-values: (zmin,zmax)"""
#     super(PeriodicBoundary3D, self).__init__()
#     self.left,  self.right = xlims
#     self.bottom,self.top   = ylims
#     self.back,  self.front = zlims
#     self.xspan = self.right-self.left
#     self.yspan = self.top-self.bottom
#     self.zspan = self.front-self.back
#   # Left boundary is "target domain"
#   def inside(self, x, on_boundary):
#     # return True if on left or bottom boundary AND NOT on one of the two corners 
#     #  (self.left, self.top) and (self.right, self.bottom)
#     return bool((fem.near(x[0], self.left) or fem.near(x[1], self.bottom)) and 
#                 (not ((fem.near(x[0], self.left) and fem.near(x[1], self.top)) or 
#                       (fem.near(x[0], self.right) and fem.near(x[1], self.bottom)))) and on_boundary)
#   def map(self, x, y):
#     if fem.near(x[0], self.right) and fem.near(x[1], self.top):
#       y[0] = x[0] - self.xspan
#       y[1] = x[1] - self.yspan
#     elif fem.near(x[0], self.right):
#       y[0] = x[0] - self.xspan
#       y[1] = x[1]
#     else:   # fem.near(x[1], self.top)
#       y[0] = x[0]
#       y[1] = x[1] - self.yspan

_HomogFickianConditions_props_schema_yaml="""#HomogFickianConditions
elementorder: {type: integer}
dirichlet:
  anyOf:
    - {type: array}
    - {type: object}
boundaries: {type: array}
"""
HomogFickianConditions_props_schema=yaml_manager.readstring(_HomogFickianConditions_props_schema_yaml)
HomogFickianConditions_schema=simrequest.update_schema_props(simrequest.EmptyConditions_schema,
                                                    HomogFickianConditions_props_schema,
                                                    ['elementorder','boundaries'])

class HomogFickianSimulator(simrequest.SimulationRequest):
  """Simulator for Homogenized Fickian Diffusion
  
  Isotropy of the input diffusion constant is assumed."""
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,HomogFickianConditions_schema)
  
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)
    
    #Mesh calculations
    spatial_dims=self.meshinfo.mesh.geometry().dim()
    meshcoords=self.meshinfo.mesh.coordinates()
    lowerlims=tuple([np.amin(meshcoords[:,i]) for i in range(spatial_dims)])
    upperlims=tuple([np.amax(meshcoords[:,i]) for i in range(spatial_dims)])
    pairedlims=list(zip(lowerlims,upperlims))

    if hasattr(conditions, 'dirichlet'):
      #The dirichlet conditions are a substitute for periodic boundary conditions
      using_periodic=False
    else:
      #Use periodic boundary conditions
      using_periodic=True
      self.bcs=[]
      if spatial_dims==2:
        pbc = PeriodicBoundary2D(*pairedlims)
      elif spatial_dims==3:
        raise NotImplementedError("Sorry, periodic 3D BCs aren't working yet.")
        #pbc = PeriodicBoundary3D(*pairedlims)
      else:
        raise NotImplementedError("You asked for a simulation on a mesh with %d spatial dimensions."%spatial_dims)

    #Function Spaces and Functions
    if using_periodic:
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder, constrained_domain=pbc)
    else:
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder)
    self.scalar_V = fem.FunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder)
    #Trial and Test Functions
    chi=fem.TrialFunction(self.V)
    v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V,name='chi')

    if hasattr(conditions,'dirichlet'):
      #Use Dirichlet boundary conditions instead of truly periodic boundary conditions
      if isinstance(conditions.dirichlet,dict):
        self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]
      else:
        #This is a workaround for the meshes that don't have meshfunctions
        val=conditions.dirichlet
        def boundary(x, on_boundary):
          ans = False
          for i in range(spatial_dims):
            for j in range(2):
              ans = ans or fem.near(x[i],pairedlims[i][j],BOUNDTOL)
          ans = ans and on_boundary
          return ans
        self.bcs=[fem.DirichletBC(self.V,val,boundary)]

    #Load diffusion constant as a Function
    if hasattr(self,'loaddata'):
      self.D=fem.Function(self.scalar_V)
      self.process_load_commands()
    else:
      self.D = fem.Constant(1)
    
    #The index objects
    i=ufl.i
    j=ufl.j

    #Measure and normal for external boundaries
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Equation term dictionary
    eqnterms=equationbuilder.EquationTermDict()

    #Bilinear boundary terms
    for psurf in conditions.boundaries:
      termname="bilinear_boundary_%d"%psurf
      form=self.D*self.n[i]*chi[j].dx(i)*v[j]*self.ds(psurf)
      eqnterms.add(termname,form,bilinear=True)

    #Bilinear body terms
    termname="bilinear_body"
    form=-self.D*chi[j].dx(i)*v[j].dx(i)*self.dx
    eqnterms.add(termname,form,bilinear=True)

    #Linear boundary terms
    for psurf in conditions.boundaries:
      termname="linear_boundary_%d"%psurf
      form=self.D*self.n[i]*v[i]*self.ds(psurf)
      eqnterms.add(termname,form,bilinear=False)

    #Linear body terms
    termname="linear_body"
    form=-self.D*v[i].dx(i)*self.dx
    eqnterms.add(termname,form,bilinear=False)

    #FEniCS Problem and Solver
    a=eqnterms.sumterms(bilinear=True)
    L=eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(a,L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

    #Get solver parameters
    self.set_solver_parameters()

    #Solve
    if not getattr(self,'skipsolve',False):
      self.solver.solve()

  def macroscale_diffusion(self,respath="D_macro",attrpath="soln",volpath="volume"):
    """Perform the integral to obtain the homogenized diffusion constant
    
    Isotropy of the input D is assumed, but the output D may be anisotropic or even non-diagonal.

    Arguments:

      - respath = optional, attribute path for storing the result
      - attrpath = optional, attribute path storing the solution (the result for chi)
      - volpath = optional, attribute path storing the unit cell volume.
          
    Required attributes (other than those from simulator_general):
    
      - the attribute given by attrpath
      - dx = FEniCS Measure for cells
      - D = FEniCS Function with the diffusion constant
    
    New attribute created/overwitten.
    No return value.
    No output files."""
    d=self.meshinfo.mesh.geometry().dim()
    volume=self.get_nested(volpath)
    kdelta = lambda i,j: 1 if i==j else 0 #Kronecker delta
    soln=self.get_nested(attrpath)
    gradchi=fem.grad(soln)
    matr=[]
    for ii in range(d):
      row=[]
      for jj in range(d):
        term1=kdelta(ii,jj)*fem.assemble(self.D*self.dx)
        term2=fem.assemble(self.D*gradchi[jj,ii]*self.dx)
        val=(term1-term2)/volume
        row.append(float(val))
      matr.append(row)
    self.set_nested(respath,matr)

#Register for loading from yaml
yaml_manager.register_classes([HomogFickianSimulator])
