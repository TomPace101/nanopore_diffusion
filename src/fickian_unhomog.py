#Solve the unhomogenized Fickian diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import os
import os.path as osp

#Site packages
import fenics as fem

#Local
from folderstructure import *
import solver_general
import useful

class UnhomogFickianConditions(solver_general.GenericConditions):
  """Condition defnitions for use with UnhomogFickianSolver
  Attributes:
    D_bulk = bulk diffusion constant"""
  __slots__=['D_bulk']

class UnhomogFickianSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized Fickian Diffusion
  Additional attributes not inherited from GenericSolver:
    conditions = instance of UnhomogFickianConditions
    V = FEniCS FunctionSpace on the mesh
    V_vec = FEniCS VectorFunctionSpace on the mesh
    bcs = FEniCS BCParameters
    ds = FEniCS Measure for surface boundary conditions
    c = FEniCS TrialFunction on V
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams):
    """Initialize the model.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""

    #Mesh setup, output init
    super().__init__(modelparams,meshparams)

    #Get conditions
    self.conditions=UnhomogFickianConditions(**modelparams.conditions)

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.mesh, "CG", self.conditions.elementorder)

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.surfaces,psurf) for psurf,val in self.conditions.bcdict.items()]

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

    #Define variational problem
    self.c=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    self.a=fem.dot(fem.grad(self.c),fem.grad(self.v))*fem.dx
    self.L=fem.Constant(0)*self.v*self.ds

  def solve(self):
    "Do the step of solving this equation"
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bcs)
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
    if modelparams.equation in solverclasses.keys():
      meshparams=allmeshes[modelparams.meshname]
      solver=solverclasses[modelparams.equation].complete(modelparams,meshparams)
