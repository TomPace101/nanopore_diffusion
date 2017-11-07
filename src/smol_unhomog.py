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

class LPBConditions(solver_general.GeneralConditions):
  """Condition defnitions for use with LPBSolver
  Attributes:
      debye_length = Debye length"""
  __slots__=['debye_length']

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
    self.lambda_D = self.conditions.debye_length
    
    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"

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

class SUConditions(solver_general.GeneralConditions):
  """Condition defnitions for use with SUSolver
  Attributes:
    D_bulk = bulk diffusion constant
    q = electric charge of ion
    beta = 1/kBT for the temperature under consideration, in units compatible with q times the potential
    potential = dictionary defining ModelParameters for electric potential
    trans_bclist = Dirichlet boundary conditions after Slotboom transformation
  Note also that the attribute bclist (inherited), contains Dirichlet conditions on c, rather than cbar.
    That is, the code will do the Slotboom transformation on the Dirichlet boundary conditions."""
  __slots__=['D_bulk','q','beta','potential','trans_bclist']
  def transform_bcs(self,potentialparams):
    """Apply Slotboom transformation to Dirichelt boundary conditions.
    This function requires that the potential and concentration have Dirichlet boundary conditions on the same surfaces.
    Arguments:
      potentialparams = solver_general.ModelParameters for the electric potential field
    No return value.
    trans_bclist attribute is updated"""
    ##TODO
    pass
    

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
    self.V = fem.FunctionSpace(self.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.mesh, "CG", self.conditions.elementorder)

    #Set up electric potential field
    potentialparams_dict=self.conditions.potential
    for key in ['modelname','meshname','meshparamsfile','basename']:
      potentialparams_dict[key]=getattr(modelparams,key)
    potentialparams=solver_general.ModelParameters(**potentialparams_dict)
    potsolv=potentialsolverclasses[potentialparams.equation].complete(potentialparams,meshparams,writeinfo=False)
    self.info['potential']=potsolv.info

    #Dirichlet boundary conditions
    self.conditions.transform_bcs(potentialparams) #apply Slotboom transformation
    self.bcs=[fem.DirichletBC(self.V,val,self.surfaces,psurf) for psurf,val in self.conditions.trans_bclist]

    #Neumann boundary conditions
    #they are all zero in this case
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.surfaces) ##TODO: specify which surfaces are Neumann?

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
