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
import unitsystem as UN
import common
import solver_general


class SpeciesInfo(common.ParameterSet):
  """Information on each species
  Attributes:
    symbol: [list of chemical symbols as strings],
    z: [list of ionic charges as numbers],
    initconc: [list of initial concentrations as numbers],
    D: [list of diffusion constants]
    N: number of species (calculated)
  Each list must have 1 entry per diffusing chemical species"""
  __slots__=['symbol','z','initconc','D','N']
  def __init__(self,**kwargs):
    #Initialization from base class
    super().__init__(**kwargs)
    #Check number of species
    nspec_all=[len(l) for l in self.conditions.species_info.values()]
    assert min(nspec_all)==max(nspec_all), "Inconsistent number of species: %s"%str(nspec_all)
    self.N=nspec_all[0]

class ReactionInfo(common.ParameterSet):
  """Information on each reaction
  Attributes:
    constants: [list of reaction rate constants],
    functions: [list of reaction rate function ....]}
    stoichio: [list of stoichiometric coefficients lists, negative for reactants, positive for products]
      each entry for stoichiometric coefficients is itself a list (one such list for each reaction), with one number for each species in the list for each reaction
      Thus, the total number of stoichiometric coefficients is the product of the number of reactions and the number of species.
  Each list must have 1 entry per uni-directional reaction (bidirectional reactions are considered as 2 uni-directional reactions each)"""
  __slots__=['constants','functions','stoichio','N']
  def __init__(self,**kwargs):
    #Initialization from base class
    super().__init__(**kwargs)
    #Check number of reactions
    nreac_all=[len(l) for l in reaction_info_dict.values()]
    assert min(nreac_all)==max(nreac_all), "Inconsistent number of reactions: %s"%str(nreac_all)    
    self.N=nreac_all[0]

class TDPNPUConditions(solver_general.GenericConditions):
  """Condition defnitions for use with TDPNPUSolver
  Attributes:
    temperature = the temperature under consideration, as a number
    species_info = dictionary defining a SpeciesInfo object
    reaction_info = dictionary defining a ReactionInfo object
    bcdict = dictionary of Dirichlet boundary conditions: {physical ...##TODO}
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
    conditions = instance of TDPNPUConditions
    species = instance of SpeciesInfo
    reactions = instance of ReactionInfo
    Nvars = number of field variables to solve for

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
    self.species=SpeciesInfo(**self.conditions.species_info)
    self.reaction=ReactionInfo(**self.conditions.reaction_info)

    ##TODO
    non_species_vars=['Phi']
    varlist=self.species.symbol+non_species_vars
    self.Nvars=len(varlist)
    
    #Function space(s)
    ##TODO

    #Measure for external boundaries
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.surfaces)

    #Dirichlet boundary conditions
    ##TODO

    #Neumann boundary conditions
    ##TODO

    #Define variational problem
    ##TODO

  def solve(self):
    "Do the time steps"
    
    ##TODO
    pass

solverclasses={'tdpnp_unhomog':TDPNPUSolver}
