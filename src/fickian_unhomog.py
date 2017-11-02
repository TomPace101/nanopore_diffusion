#Solve the unhomogenized Fickian diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import os
import os.path as osp

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
  def to_bclist(self,fs,surfaces):
    """Output list of FEniCS DirichletBC objects based on given parameters
    Arguments:
      fs = FEniCS FunctionSpace to use for the boundary conditions
      surfaces = FEniCS MeshFunction to use for the surfaces
    Returns:
      bcs = list of DirichletBC objects"""
    bcs=[]
    dpairs=[(self.basesurf,self.baseval), (self.topsurf,self.topval)] #Physical surface and Dirichlet value pairs
    for psurf,val in dpairs:
      bcs.append(DirichletBC(fs,val,surfaces,psurf))
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
  def __init__(self,modelparams,meshparams):
    """Initialize the model, and optionally solve and generate output.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
    
    #Mesh setup
    super().__init__(modelparams,meshparams)
    
    #Function space for scalars and vectors
    self.V = FunctionSpace(self.mesh,'CG',1) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = VectorFunctionSpace(self.mesh, "CG", 1)

    #Dirichlet boundary conditions
    self.bcs=BCParameters(**self.modelparams.conditions).to_bclist(self.V, self.surfaces)

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

    #Define variational problem
    self.c=TrialFunction(self.V)
    self.v=TestFunction(self.V)
    self.a=dot(grad(self.c),grad(self.v))*dx
    self.L=Constant(0)*self.v*self.ds
    
  def solve(self):
    "Do the step of solving this equation"
    self.soln=Function(self.V)
    solve(self.a==self.L, self.soln, self.bcs)
    return

solverclasses={'fickian_unhomog':UnhomogFickianSolver}

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Solve the unhomogenized fickian diffusion equation with fenics')
  parser.add_argument('model_params_file', help='filename (not complete path) containing ModelParameters definitions')
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.model_params_file), "Model parameter definition file does not exist: %s"%(cmdline.model_params_file)

  #Get all models to solve, and all their meshes
  allmodels,modelfiles,allmeshes,meshfiles=solver_general.GetAllModelsAndMeshes([cmdline.model_params_file])
  
  #Run each requested analysis
  for modelparams in allmodels.values():
    #Only do analyses with equations supported by this module
    if modelparams.equation in solverclasses.keys() 
      meshparams=allmeshes[modelparams.meshname]
      solver=solverclasses[modelparams.equation].complete(modelparams,meshparams)

