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
      modelparams = solver_run.ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
      
      ##TODO
      pass

  def solve(self):
    "Do the step of solving this equation"
    
    ##TODO
    pass

solverclasses={'pnp_unhomog':PNPUSolver}