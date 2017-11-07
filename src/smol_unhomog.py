#Solve the unhomogenized Smoluchowski diffusion problem
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

class LPBConditions(useful.ParameterSet):
  """Condition defnitions for use with LPBSolver
  Attributes:
    bclist = list of Dirichlet boundary conditions pairs:
      [(physical surface number, solution value), ...]"""
  __slots__=['bclist']

class LPBSolver(solver_general.GenericSolver):
  """Solver for linearized Poisson-Boltzmann equation.
  Additional attributes not inherited from GenericSolver:
    conditions = instance of LPBConditions
    lambda_D = Debye length
    V = FEniCS FunctionSpace on the mesh
    bcs = list of FEniCS DirichletBC instances
    ds = FEniCS Measure for surface boundary conditions
    phi = FEniCS TrialFunction on V
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams):
    """Initialize the model, and optionally solve and generate output.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance
      complete = boolean, True to solve and generate output"""
    #Mesh setup
    super().__init__(modelparams,meshparams)
    
    #Get conditions
    self.conditions=LPBConditions(**modelparams.conditions)
    
    #Properties of problem domain
    self.lambda_D = modelparams.properties['debye_length']
    
    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.mesh,'CG',modelparams.elementorder) #CG="continuous galerkin", ie "Lagrange"

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.surfaces,psurf) for psurf,val in self.conditions.bclist]

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

    #Define variational problem
    self.phi=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    self.a=((1/self.lambda_D**2)*self.phi*self.v + fem.dot(fem.grad(self.phi),fem.grad(self.v)))*fem.dx
    self.L=fem.Constant(0)*self.v*self.ds
    
  def solve(self):
    "Do the step of solving this equation"
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bcs)
    return

#Lookup of electric potential solvers by name
potentialsolverclasses={'linear_pb':LPBSolver}

class SUConditions(useful.ParameterSet):
  """Condition defnitions for use with SUSolver
  Attributes:
    bclist = list of Dirichlet boundary conditions pairs:
      [(physical surface number, solution value), ...]
    potential = dictionary defining ModelParameters for electric potential""" ##TODO: is solution value in bclist c or cbar?
  __slots__=['bclist','potential']

class SUSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized Smoluchowski Diffusion
  Additional attributes not inherited from GenericSolver:
    conditions = instance of SUConditions
    V = FEniCS FunctionSpace on the mesh
    V_vec = FEniCS VectorFunctionSpace on the mesh
    bcs = FEniCS BCParameters
    ds = FEniCS Measure for surface boundary conditions
    c = FEniCS TrialFunction on V
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams,meshparams):
    """Initialize the model, and optionally solve and generate output.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
    
    #Mesh setup, output init
    super().__init__(modelparams,meshparams)
    
    #Get conditions
    self.conditions=SUConditions(**modelparams.conditions)

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.mesh,'CG',modelparams.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.mesh, "CG", modelparams.elementorder)

    #Dirichlet boundary conditions
    self.bcs=[fem.DirichletBC(self.V,val,self.surfaces,psurf) for psurf,val in self.conditions.bclist]

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

    #Set up electric potential field
    potentialparams_dict=self.conditions.potential
    for key in ['modelname','meshname','meshparamsfile','basename']:
      potentialparams_dict[key]=getattr(modelparams,key)
    potentialparams=solver_general.ModelParameters(**potentialparams_dict)
    potsolv=potentialsolverclasses[potentialparams.equation].complete(potentialparams,meshparams,writeinfo=False)
    self.info['potential']=potsolv.info

    #Define variational problem
    self.c=fem.TrialFunction(self.V)
    self.v=fem.TestFunction(self.V)
    ##TODO: input smoluchowski weak form here
    self.a=fem.dot(fem.grad(self.c),fem.grad(self.v))*fem.dx ##TODO
    self.L=fem.Constant(0)*self.v*self.ds ##TODO
    
  def solve(self):
    "Do the step of solving this equation"
    ##TODO: slotboom!
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bcs)
    return

solverclasses={'smol_unhomog':SUSolver}

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
    #Only do analyses with equations supported by this module
    if modelparams.equation in solverclasses.keys():
      meshparams=allmeshes[modelparams.meshname]
      ##TODO: uncomment when ready
      ##solver=solverclasses[modelparams.equation].complete(modelparams,meshparams)
      solver=solverclasses[modelparams.equation](modelparams,meshparams)
