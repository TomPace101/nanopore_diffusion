"""Solve the homogenized Fickian diffusion problem in 3D
and extract the data needed for post-processing efforts"""

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

BOUNDTOL=1e-6

##TODO
#Modify the class below to be 3D
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

_HomogFickian3DConditions_props_schema_yaml="""#HomogFickian3DConditions
elementorder: {type: integer}
dirichlet: {type: object}
boundaries: {type: array}
"""
HomogFickian3DConditions_props_schema=yaml_manager.readstring(_HomogFickian3DConditions_props_schema_yaml)
HomogFickian3DConditions_schema=simrequest.update_schema_props(simrequest.EmptyConditions_schema,
                                                    HomogFickian3DConditions_props_schema,
                                                    ['elementorder','boundaries'])

class HomogFickian3DSimulator(simrequest.SimulationRequest):
  """Simulator for 3D Homogenized Fickian Diffusion
  
  Isotropy of the input diffusion constant is assumed."""
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,HomogFickian3DConditions_schema)
  
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Function Spaces and Functions
    #Function spaces
    # self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder, constrained_domain=pbc)
    self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder) #Temporary: use Dirichlet conditions instead
    self.scalar_V = fem.FunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder)
    #Trial and Test Functions
    chi=fem.TrialFunction(self.V)
    v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V,name='chi')

    #Periodic boundary condition
    spatial_dims=self.meshinfo.mesh.geometry().dim()
    Npts=self.scalar_V.dim()
    # xkeys=[k for k in self.meshinfo.metadata.keys() if k[0].upper()=='X']
    # ykeys=[k for k in self.meshinfo.metadata.keys() if k[0].upper()=='Y']
    # zkeys=[k for k in self.meshinfo.metadata.keys() if k[0].upper()=='Z']
    # xvals=[self.meshinfo.metadata[k] for k in xkeys]
    # yvals=[self.meshinfo.metadata[k] for k in ykeys]
    # zvals=[self.meshinfo.metadata[k] for k in zkeys]
    # xlims=(min(xvals),max(xvals))
    # ylims=(min(yvals),max(yvals))
    # zlims=(min(zvals),max(zvals))
    ptcoords=self.scalar_V.tabulate_dof_coordinates().reshape(Npts,spatial_dims)
    lowerlims=tuple([np.amin(ptcoords[:,i]) for i in range(spatial_dims)])
    upperlims=tuple([np.amax(ptcoords[:,i]) for i in range(spatial_dims)])
    pairedlims=list(zip(lowerlims,upperlims))
    xlims,ylims,zlims=pairedlims
    # pbc = PeriodicBoundary3D(xlims,ylims,zlims)

    #Dirichlet boundary conditions
    if isinstance(conditions.dirichlet,dict):
      self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]
    else:
      #This is a temporary workaround for the zeolite meshes that don't have meshfunctions
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
    i=fem.i
    j=fem.j

    #Measure and normal for external boundaries
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Equation term dictionary
    eqnterms=equationbuilder.EquationTermDict()

    #Bilinear boundary terms
    for psurf in conditions.boundaries:
      termname="bilinear_boundary_%d"%psurf
      ufl=self.D*self.n[i]*chi[j].dx(i)*v[j]*self.ds(psurf)
      eqnterms.add(termname,ufl,bilinear=True)

    #Bilinear body terms
    termname="bilinear_body"
    ufl=-self.D*chi[j].dx(i)*v[j].dx(i)*self.dx
    eqnterms.add(termname,ufl,bilinear=True)

    #Linear boundary terms
    for psurf in conditions.boundaries:
      termname="linear_boundary_%d"%psurf
      ufl=self.D*self.n[i]*v[i]*self.ds(psurf)
      eqnterms.add(termname,ufl,bilinear=False)

    #Linear body terms
    termname="linear_body"
    ufl=-self.D*v[i].dx(i)*self.dx
    eqnterms.add(termname,ufl,bilinear=False)

    #FEniCS Problem and Solver
    a=eqnterms.sumterms(bilinear=True)
    L=eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(a,L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

    #Solve
    self.solver.parameters['linear_solver']='gmres'
    self.solver.parameters['preconditioner']='ilu'
    self.solver.parameters['krylov_solver']['absolute_tolerance']=1e-12
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
    volume=self.get_nested(volpath)
    kdelta = lambda i,j: 1 if i==j else 0 #Kronecker delta
    soln=self.get_nested(attrpath)
    gradchi=fem.grad(soln)
    matr=[]
    for ii in range(3):
      row=[]
      for jj in range(3):
        term1=kdelta(ii,jj)*fem.assemble(self.D*self.dx)
        term2=fem.assemble(self.D*gradchi[jj,ii]*self.dx)
        val=(term1-term2)/volume
        row.append(val)
      matr.append(row)
    self.set_nested(respath,matr)

#Register for loading from yaml
yaml_manager.register_classes([HomogFickian3DSimulator])
