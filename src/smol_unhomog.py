#Solve the unhomogenized Smoluchowski diffusion problem
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

class SUSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized Smoluchowski Diffusion
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
    self.V = FunctionSpace(self.mesh,'CG',1) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = VectorFunctionSpace(self.mesh, "CG", 1)

    #Dirichlet boundary conditions
    self.bcs=BCParameters(**self.modelparams.conditions).to_bclist(self.V, self.surfaces)

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

    #Set up electric potential field
    ##TODO

    #Define variational problem
    self.c=TrialFunction(self.V)
    self.v=TestFunction(self.V)
    ##TODO: input smoluchowski weak form here
    self.a=dot(grad(self.c),grad(self.v))*dx ##TODO
    self.L=Constant(0)*self.v*self.ds ##TODO
    
    #If requested, solve and generate output
    if complete:
      self.complete()

  def solve(self):
    "Do the step of solving this equation"
    self.soln=Function(self.V)
    solve(self.a==self.L, self.soln, self.bcs)
    return

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Solve the unhomogenized Smoluchowski diffusion equation with fenics')
  parser.add_argument('model_params_file', help='filename (not complete path) containing ModelParameters definitions')
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.model_params_file), "Model parameter definition file does not exist: %s"%(cmdline.model_params_file)

  #Get all models to solve, and all their meshes
  allmodels,modelfiles,allmeshes,meshfiles=solver_general.GetAllModelsAndMeshes([cmdline.model_params_file])
  
  #Run each requested analysis
  for modelparams in allmodels.values():
    meshparams=allmeshes[modelparams.meshname]
    solver=SUSolver(modelparams,meshparams)
    solver.complete()
