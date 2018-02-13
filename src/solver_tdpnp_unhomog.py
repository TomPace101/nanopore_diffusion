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
import rxn_rate_funcs


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
    functions: [list of reaction rate function ....] (all functions must be in the reaction rates module)
    stoichio: [list of stoichiometric coefficients lists, negative for reactants, positive for products]
      each entry for stoichiometric coefficients is itself a list (one such list for each reaction), with one number for each species in the list for each reaction
      Thus, the total number of stoichiometric coefficients is the product of the number of reactions and the number of species.
    N: number of reactions (calculated)
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
    elementorder = see base class
    bcdict = dictionary of Dirichlet boundary conditions:
      {physical facet number: [solution values, None for no condition]}
    temperature = the temperature under consideration, as a number
    eps_r = relative permittivity of the medium
    species_info = dictionary defining a SpeciesInfo object
    reaction_info = dictionary defining a ReactionInfo object
    initial_potential = initial electric potential, assumed constant over space, as number
    t_end = end time for simulation (may be exceeded if not exactly divisible by timestep)
    delta_t = timestep for simulation (number of timesteps is calculated from this)
    beta = optional, calculated from temperature if not provided"""
  __slots__=['beta','temperature','eps_r','species_info','reaction_info','initial','t_end','delta_t']
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
    Nspecies = number of chemical species
    Nvars = number of field variables to solve for
    dt = timestep
    numsteps = number of steps to compute
    V = FEniCS FunctionSpace on the mesh
    bcs = FEniCS BCParameters
    ds = FEniCS Measure for facet boundary conditions
    FF = symbolic functional form, which is set equal to zero in the weak form equation
    J = symbolic Jacobian of FF

    v = FEniCS TestFunction on V"""
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
    self.reactions=ReactionInfo(**self.conditions.reaction_info)

    #List and count the degrees of freedom
    self.Nspecies=len(species.symbol)
    non_species_vars=['Phi']
    varlist=self.species.symbol+non_species_vars
    self.Nvars=len(varlist)
    
    #Calculate number of time steps
    self.dt=self.conditions.delta_t
    self.numsteps=math.ceil(t_end/self.dt)
    
    #Elements and Function space(s)
    ele = fem.FiniteElement('P',mesh.ufl_cell(),elementorder)
    mele = fem.MixedElement([ele]*Nvars)
    self.V = fem.FunctionSpace(self.mesh,mele)

    #Test and trial functions
    u = fem.Function(V)
    trialfuncs=fem.split(u)
    clist=trialfuncs[:self.Nspecies]
    Phi=trialfuncs[self.Nspecies]
    testfuncs=fem.TestFunctions(self.V)
    vlist=testfuncs[:self.Nspecies]
    vPhi=testfuncs[self.Nspecies]

    #Measure for external boundaries
    self.ds = fem.Measure("ds",domain=self.mesh,subdomain_data=self.facets)

    #Dirichlet boundary conditions
    self.bcs=[]
    for psurf,vals in bcdict.items():
      for i,value in enumerate(vals):
        if value is not None:
          self.bcs.append(fem.DirichletBC(V.sub(i),Constant(value),self.facets,psurf))

    #Neumann boundary conditions
    ##TODO

    #Initial Conditions and Guess
    guesstup=self.species.initconc+[self.conditions.initial_potential]
    u.interpolate(fem.Constant(guesstup))
    u_k=fem.interpolate(fem.Constant(guesstup),self.V)
    u_klist=split(u_k)
    c_klist=u_klist[:self.Nspecies]
    
    #Weak Form
    #Steady-state Nernst-Planck terms for each species
    ##TODO: use revised weak form
    weakforms=[]
    for i,c in enumerate(clist):
      if self.species.D[i] is not None:
        J=-self.species.D[i]*(fem.grad(c)+self.conditions.beta*self.species.z[i]*c*fem.grad(Phi))
        wkform=inner(self.J,fem.grad(vlist[i]))*fem.dx
        weakforms.append(wkform)
    #Poisson terms
    poissonL=self.inner(UN.eps_0*self.conditions.eps_r*fem.grad(Phi),grad(vPhi))*fem.dx
    terms=[self.species.z[i]*clist[i] for i in range(self.species.N)]
    poissonR=sum(terms)*vPhi*fem.dx
    #Add up to Steady-State PNP
    FF_ss=sum(weakforms)+poissonL-poissonR
    #Time-dependent terms
    tdweaks=[]
    for i,c in enumerate(clist):
      term1=c*vlist[i]*dx
      term2=c_klist[i]*vlist[i]*dx
      tdweaks.append(term1-term2)
    #Reaction terms
    rxnweaks=[]
    for i,c in enumerate(clist):
      for j in range(self.reactions.N):
        if self.reactions.stoichio[j][i] != 0:
          termconst=self.reactions.stoichio[j][i]*self.reactions.constants[j]
          rxf=getattr(rxn_rate_funcs,self.reactions.functions[j])
          term=termconst*rxf(*clist)*vlist[i]*dx
          rxnweaks.append(term)
    #Put it all together
    self.FF=sum(tdweaks)-self.dt*FF_ss-self.dt*sum(rxnweaks)
    #Take derivative
    self.J=derivative(self.FF,u)

  def solve(self):
    "Do the time steps"
    
    ##TODO: open the output files
    for k in range(self.numsteps):
      ##TODO
      pass

solverclasses={'tdpnp_unhomog':TDPNPUSolver}