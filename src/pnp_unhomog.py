#Solve the unhomogenized PNP diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import math
import os
import os.path as osp

#Site packages
import fenics as fem

#Local
from folderstructure import *
import solver_general
import useful

class PNPUConditions(solver_general.GenericConditions):
  """Condition defnitions for use with PNPUSolver
  Attributes:
    beta = 1/kBT for the temperature under consideration, in units compatible with q times the potential
    species_info = dictionary {symbol: charge}"""
  __slots__=['beta','species_info']

class PNPUSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized PNP Diffusion
  Additional attributes not inherited from GenericSolver:
    conditions = instance of PNPUConditions
    V = FEniCS FunctionSpace on the mesh
    V_vec = FEniCS VectorFunctionSpace on the mesh
    bcs = FEniCS BCParameters
    ds = FEniCS Measure for surface boundary conditions
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams):
    """Initialize the model.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
      
      ##TODO
      pass

  def solve(self):
    "Do the step of solving this equation"
    
    ##TODO
    pass

solverclasses={'pnp_unhomog':PNPUSolver}

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Solve the unhomogenized PNP diffusion equation with fenics')
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

