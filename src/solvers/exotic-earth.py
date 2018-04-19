#Solve the homogenized Fickian diffusion problem, specific to the analysis needed by exotic-earth

#Standard Library
from argparse import Namespace

#Site Packages
import numpy as np
import fenics as fem

#Local
import solver_general

class VaryingPropertyByCell(fem.Expression):
  """TODO: document"""
  def __init__(self,meshfunc,prop_map,**kwargs):
    self.meshfunc=meshfunc
    self.prop_map=prop_map
  def value_shape(self):
    return (2,2)
  def eval_cell(self, values, x, cell):
    key=self.meshfunc[cell.index]
    for i in range(2):
      for j in range(2):
        values[2*j+i]=self.prop_map[key][i,j]

##TODO: replace with fully periodic boundary condition using `near`
class PeriodicBoundary(fem.SubDomain):
  "Copied from FEniCS example"
  # Left boundary is "target domain" G
  def inside(self, x, on_boundary):
    return bool(x[0] < fem.DOLFIN_EPS and x[0] > -fem.DOLFIN_EPS and on_boundary)
  # Map right boundary (H) to left boundary (G)
  def map(self, x, y):
    y[0] = x[0] - 1.0
    y[1] = x[1]

class ExoticEarthConditions(solver_general.GenericConditions):
  """Condition defnitions for use with ExoticEarthSolver
  Attributes:
    dirichlet = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}
    isotropic_D_values = pair of diffusion constants: regions 1 and 2, respectively"""
  __slots__=['dirichlet','isotropic_D_values']

class ExoticEarthSolver(solver_general.GenericSolver):
  """Solver for Homogenized Fickian Diffusion, special case
  Additional attributes not inherited from GenericSolver:
    conditions = instance of ExoticEarthConditions
    V = FEniCS FunctionSpace on the mesh
    bclist = list of FEniCS DirichletBC instances
    ds = FEniCS Measure for facet boundary conditions
    chi = FEniCS TrialFunction on V
    v = FEniCS TestFunction on V
    a = bilinear form in variational problem
    L = linear form in variational problem"""
  def __init__(self,modelparams):
    """Initialize the model.
    Arguments:
      modelparams = solver_run.ModelParameters instance"""

    #Load parameters, init output, mesh setup
    super(ExoticEarthSolver, self).__init__(modelparams)
    self.loadmesh()

    #Get conditions
    self.conditions=UnhomogFickianConditions(**modelparams.conditions)
    
    #Isotropic diffusion matrices of the two regions
    D1,D2=self.conditions.isotropic_D_values
    D1arr=np.array([[D1,0],[0,D1]])
    D2arr=np.array([[D2,0],[0,D2]])
    D_prop_map={1:D1arr,2:D2arr}

    #Geometric properties based on parametric locations
    parmlocs=self.parametric_locations
    Ly = paramlocs.Y3-paramlocs.Y1
    Lx = paramlocs.X2-paramlocs.X1
    alpha1=(paramlocs.Y3-paramlocs.Y2)/Ly
    alpha2=(paramlocs.Y2-paramlocs.Y1)/Ly
    
    #Theoretical result
    Dxx_theo=alpha1*D1+alpha2*D2
    Dyy_theo=D1*D2/(alpha1*D2+alpha2*D1)
    Dtheo_arr=np.array([[Dxx_theo,0],[0,Dyy_theo]])

    #Varying D
    Darr=VaryingPropertyByCell(cells,D_prop_map,degree=0)

    # Create periodic boundary condition
    pbc = PeriodicBoundary()

    #Function space and functions for weak form
    self.V = fem.VectorFunctionSpace(mesh, 'P', elementorder, constrained_domain=pbc)
    self.chi=fem.TrialFunction(V)
    self.v=fem.TestFunction(V)


    #Dirichlet Boundary Conditions
    self.bclist=[]
    dirichlet=getattr(self.conditions,'dirichlet',{})
    for psurf,tup in dirichlet.items():
      vec=fem.Constant(tup)
      self.bclist.append(fem.DirichletBC(self.V,vec,self.facets,psurf))

    #Weak Form
    gradmat=fem.dot(fem.grad(self.v).T,fem.grad(self.chi))
    A=fem.inner(Darr,gradmat)*fem.dx
    L=fem.inner(Darr,fem.grad(self.v).T)*fem.dx

  def solve(self):
    "Do the step of solving this equation"
    self.soln=fem.Function(self.V)
    fem.solve(self.a==self.L, self.soln, self.bclist)
    return

 def calc_Dmacro(self):
    """Do the integral for the effective diffusion constant"""
    #TODO: documentation
    d3x = fem.Measure('cell',domain=self.mesh)
    volume=fem.assemble(fem.Constant(1)*d3x)
    int_mat=fem.dot(Darr,fem.grad(self.soln).T)
    matr=[]
    for i in range(2):
      row=[]
      for j in range(2):
        val=(fem.assemble(Darr[i,j]*d3x)-fem.assemble(int_mat[i,j]*d3x))/volume
        row.append(val)
      matr.append(row)
    self.Dmacro=np.array(matr)

solverclasses={'exotic-earth':ExoticEearthSolver}
