#Solve the unhomogenized Fickian diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import os
import os.path as osp
import sys

#Site packages
from fenics import *

#Local
from folderstructure import *
import solver_general
import useful

class BCParameters(useful.ParameterSet):
  """Boundary condition definition
  Attributes:
    topsurf = physical surface number for top surface
    basesurf = physical surface number for base surface
    topval = concentration value at top surface
    baseval = concentration value at base surface"""
  __slots__=['topsurf','basesurf','topval','baseval']
  def to_bclist():
    "output list of FEniCS DirichletBC objects based on given parameters"
    bcs=[]
    dpairs=[(self.basesurf,self.baseval), (self.topsurf,self.topval)] #Physical surface and Dirichlet value pairs
    for psurf,val in dpairs:
      bcs.append(DirichletBC(V,val,surfaces,psurf))
    return bcs

class UnhomogFickianSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized Fickian Diffusion
  Additional attributes not inherited from GenericSolver:
    V = FEniCS FunctionSpace on the mesh
    V_vec = FEniCS VectorFunctionSpace on the mesh
    bcs = FEniCS BCParameters
    ds = FEniCS Measure for surface boundary conditions
    c = FEniCS TrialFunction on V
    v = FEniCS TestFunctoin on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams,complete=False):
    """Initialize the model, and optionally solve and generate output.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance
      complete = boolean, True to solve and generate output"""
    
    #Mesh setup
    super().__init__(modelparams,meshparams)
    
    #Function space for scalars and vectors
    self.V = FunctionSpace(mesh,'CG',1) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = VectorFunctionSpace(mesh, "CG", 1)

    #Dirichlet boundary conditions
    self.bcs=BCParameters(**self.modelparams.boundaryconditions).to_bclist()
    
    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

    #Define variational problem
    self.c=TrialFunction(V)
    self.v=TestFunction(V)
    self.a=dot(grad(c),grad(v))*dx
    self.L=Constant(0)*v*ds
    
    #If requested, solve and generate output
    if complete:
      self.complete()

  def solve(self):
    "Do the step of solving this equation"
    self.soln=Function(V)
    solve(a==L,self.soln,bcs)
    return

#Support command-line arguments
if __name__ == '__main__':
  ##TODO: this needs to be updated
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Solve the unhomogenized fickian diffusion equation with fenics')
  parser.add_argument('bc_params_yaml', help='path to boundary conditions parameter yaml file')
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.bc_params_yaml), "Boundary conditions parameter definition file does not exist: %s"%cmdline.bc_params_yaml

  #Read in the yaml file
  solruns=useful.readyaml_multidoc(cmdline.bc_params_yaml)
  
  #Run each requested analysis
  for run in solruns:
    params=argparse.Namespace(**run)
    SolveMesh(params)

