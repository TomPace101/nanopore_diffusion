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
import folderstructure as FS
import unitsystem as UN
import solver_general


class TDPNPUConditions(solver_general.GenericConditions):
  """Condition defnitions for use with TDPNPUSolver
  Attributes:
    temperature = the temperature under consideration, as a number
    species_info = dictionary
      {symbol: [list of chemical symbols as strings],
       z: [list of ionic charges as numbers],
       initconc: [list of initial concentrations as numbers],
       D: [list of diffusion constants]}
      each list must have 1 entry per diffusing chemical species
    reaction_info = dictionary
      {constants: [list of reaction rate constants],
       functions: [list of reaction rate function ....]}
       stoichio: [list of stoichiometric coefficients lists, negative for reactants, positive for products]
        each entry for stoichiometric coefficients is itself a list (one such list for each reaction), with one number for each species in the list for each reaction
        Thus, the total number of stoichiometric coefficients is the product of the number of reactions and the number of species.
      each list must have 1 entry per uni-directional reaction (bidirectional reactions are considered as 2 uni-directional reactions each)
    beta = optional, calculated from temperature if not provided"""
  __slots__=['beta','temperature','species_info','reaction_info']
  def __init__(self,**kwargs):
    #Initialization from base class
    super().__init__(**kwargs)
    #If beta not provided, calculate from temperature
    if not hasattr(self,'beta'):
      self.beta = 1.0/(self.temperature*UN.kB)

class TDPNPUSolver(solver_general.GenericSolver):
  """Solver for Unhomogenized Time-Domain Poisson-Nernst-Planck Diffusion
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
      
    #Load parameters, init output, mesh setup
    super().__init__(modelparams,meshparams)
    self.loadmesh()

    #Get conditions
    self.conditions=TDPNPUConditions(**modelparams.conditions)

    #Function space for scalars and vectors
    self.V = fem.FunctionSpace(self.mesh,'CG',self.conditions.elementorder) #CG="continuous galerkin", ie "Lagrange"
    self.V_vec = fem.VectorFunctionSpace(self.mesh, "CG", self.conditions.elementorder)

    #Measure for external boundaries
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.surfaces)

    #Dirichlet boundary conditions
    self.conditions.transform_bcs(self.potsolv.conditions.bcdict,self.beta_q) #apply Slotboom transformation
    self.bcs=[fem.DirichletBC(self.V,val,self.surfaces,psurf) for psurf,val in self.conditions.trans_bcdict.items()]

    #Neumann boundary conditions
    ##TODO

    #Define variational problem
    self.d3x = fem.Measure('cell',domain=self.mesh)
    ##TODO

  def solve(self):
    "Do the time steps"
    
    ##TODO
    pass

solverclasses={'tdpnp_unhomog':TDPNPUSolver}
