"""Solve the homogenized Smoluchowski diffusion problem (corrector problem)
and extract the data needed for post-processing efforts"""

#Standard library
from argparse import Namespace

#Site packages
import numpy as np
import fenics as fem
import ufl

#This package
from ..requesthandler import yaml_manager
from . import simrequest
from . import equationbuilder
from . import periodic_boundaries

BOUNDTOL=1e-6

_HomogSmolConditions_props_schema_yaml="""#HomogSmolConditions
elementorder: {type: integer}
dirichlet:
  anyOf:
    - {type: array}
    - {type: object}
beta: {type: number}
boundaries: {type: array}
"""
HomogSmolConditions_props_schema=yaml_manager.readstring(_HomogSmolConditions_props_schema_yaml)
HomogSmolConditions_schema=simrequest.update_schema_props(simrequest.EmptyConditions_schema,
                                                    HomogSmolConditions_props_schema,
                                                    ['elementorder','boundaries','beta'])

class HomogSmolSimulator(simrequest.SimulationRequest):
  """Simulator for Homogenized Smoluchowski Diffusion
  
  Isotropy of the input diffusion constant is assumed."""
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,HomogSmolConditions_schema)
  
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)
    
    #Mesh calculations
    spatial_dims=self.meshinfo.mesh.geometry().dim()
    meshcoords=self.meshinfo.mesh.coordinates()
    lowerlims=tuple([np.amin(meshcoords[:,i]) for i in range(spatial_dims)])
    upperlims=tuple([np.amax(meshcoords[:,i]) for i in range(spatial_dims)])
    pairedlims=list(zip(lowerlims,upperlims))

    if hasattr(conditions, 'dirichlet'):
      #The dirichlet conditions are a substitute for periodic boundary conditions
      using_periodic=False
    else:
      #Use periodic boundary conditions
      using_periodic=True
      self.bcs=[]
      if spatial_dims==2:
        pbc = periodic_boundaries.PeriodicBoundary2D(*pairedlims)
      elif spatial_dims==3:
        raise NotImplementedError("Sorry, periodic 3D BCs aren't working yet.")
        #pbc = periodic_boundaries.PeriodicBoundary3D(*pairedlims)
      else:
        raise NotImplementedError("You asked for a simulation on a mesh with %d spatial dimensions."%spatial_dims)

    #Function Spaces and Functions
    if using_periodic:
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder, constrained_domain=pbc)
    else:
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder)
    self.scalar_V = fem.FunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder)
    #Trial and Test Functions
    chi=fem.TrialFunction(self.V)
    v=fem.TestFunction(self.V)
    #Solution function
    self.soln=fem.Function(self.V,name='chi')

    if hasattr(conditions,'dirichlet'):
      #Use Dirichlet boundary conditions instead of truly periodic boundary conditions
      if isinstance(conditions.dirichlet,dict):
        self.bcs=[fem.DirichletBC(self.V,val,self.meshinfo.facets,psurf) for psurf,val in conditions.dirichlet.items()]
      else:
        #This is a workaround for the meshes that don't have meshfunctions
        val=conditions.dirichlet
        def boundary(x, on_boundary):
          ans = False
          for i in range(spatial_dims):
            for j in range(2):
              ans = ans or fem.near(x[i],pairedlims[i][j],BOUNDTOL)
          ans = ans and on_boundary
          return ans
        self.bcs=[fem.DirichletBC(self.V,val,boundary)]

    #Load diffusion constant and potential as Functions
    if hasattr(self,'loaddata'):
      self.D=fem.Function(self.scalar_V)
      self.potential=fem.Function(self.scalar_V)
      self.process_load_commands()
    
    #Apply Slotboom transformation
    Dbar_expr=self.D*fem.exp(-conditions.beta*self.potential)
    self.Dbar=fem.project(Dbar_expr,self.scalar_V,solver_type="cg",preconditioner_type="amg") #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)
    self.Dbar.rename('Dbar','transformed diffusion coefficient')
    
    #The index objects
    i=ufl.i
    j=ufl.j

    #Measure and normal for external boundaries
    self.ds = fem.Measure('exterior_facet',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    self.n = fem.FacetNormal(self.meshinfo.mesh)
    self.dx = fem.Measure('cell',domain=self.meshinfo.mesh)

    #Equation term dictionary
    eqnterms=equationbuilder.EquationTermDict()

    #Bilinear boundary terms
    for psurf in conditions.boundaries:
      termname="bilinear_boundary_%d"%psurf
      form=self.Dbar*self.n[i]*chi[j].dx(i)*v[j]*self.ds(psurf)
      eqnterms.add(termname,form,bilinear=True)

    #Bilinear body terms
    termname="bilinear_body"
    form=-self.Dbar*chi[j].dx(i)*v[j].dx(i)*self.dx
    eqnterms.add(termname,form,bilinear=True)

    #Linear boundary terms
    for psurf in conditions.boundaries:
      termname="linear_boundary_%d"%psurf
      form=self.Dbar*self.n[i]*v[i]*self.ds(psurf)
      eqnterms.add(termname,form,bilinear=False)

    #Linear body terms
    termname="linear_body"
    form=-self.Dbar*v[i].dx(i)*self.dx
    eqnterms.add(termname,form,bilinear=False)

    #FEniCS Problem and Solver
    self.a=eqnterms.sumterms(bilinear=True)
    self.L=eqnterms.sumterms(bilinear=False)
    problem=fem.LinearVariationalProblem(self.a,self.L,self.soln,bcs=self.bcs)
    self.solver=fem.LinearVariationalSolver(problem)

    #Get solver parameters
    self.set_solver_parameters()

    #Solve
    if not getattr(self,'skipsolve',False):
      self.solver.solve()

  def macroscale_diffusion(self,respath="D_macro",attrpath="soln",volpath="volume"):
    """Perform the integral to obtain the homogenized diffusion constant
    
    No Slotboom transformation is applied to the result.
    (That is, the result is technically a Dbar, not just a D.)

    Isotropy of the input Dbar is assumed, but the output Dbar may be anisotropic or even non-diagonal.

    Arguments:

      - respath = optional, attribute path for storing the result
      - attrpath = optional, attribute path storing the solution (the result for chi)
      - volpath = optional, attribute path storing the unit cell volume.
          
    Required attributes (other than those from simulator_general):
    
      - the attribute given by attrpath
      - dx = FEniCS Measure for cells
      - Dbar = FEniCS Function with the resulting diffusion constant
    
    New attribute created/overwitten.
    No return value.
    No output files."""
    d=self.meshinfo.mesh.geometry().dim()
    volume=self.get_nested(volpath)
    kdelta = lambda i,j: 1 if i==j else 0 #Kronecker delta
    soln=self.get_nested(attrpath)
    gradchi=fem.grad(soln)
    matr=[]
    for ii in range(d):
      row=[]
      for jj in range(d):
        term1=kdelta(ii,jj)*fem.assemble(self.Dbar*self.dx)
        term2=fem.assemble(self.Dbar*gradchi[jj,ii]*self.dx)
        val=(term1-term2)/volume
        row.append(float(val))
      matr.append(row)
    self.set_nested(respath,matr)

#Register for loading from yaml
yaml_manager.register_classes([HomogSmolSimulator])
