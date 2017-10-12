#Solve the unhomogenized Fickian diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import os
import os.path as osp
import sys

#Site packages
from fenics import *

#Local
from folderstructure import *
import solver_general
import useful

class BCParameters(useful.ParameterSet):
  """Boundary condition definition
  Attributes:
    topsurf = physical surface number for top surface
    basesurf = physical surface number for base surface
    topval = concentration value at top surface
    baseval = concentration value at base surface"""
  __slots__=['topsurf','basesurf','topval','baseval']
  def to_bclist():
    "output list of FEniCS DirichletBC objects based on given parameters"
    bcs=[]
    dpairs=[(self.basesurf,self.baseval), (self.topsurf,self.topval)] #Physical surface and Dirichlet value pairs
    for psurf,val in dpairs:
      bcs.append(DirichletBC(V,val,surfaces,psurf))
    return bcs

#TODO: there are probably parts of this that should be refactored into functions in a more general file, once one exists
def SolveMesh(modelparams, meshparams):
  """Solve the unhomogenized Fickian diffusion equation on the indicated mesh.
  Arguments:
    modelparams = solver_general.ModelParameters instance
      dataextraction:
        fluxsurf: physical surface number for flux measurement
        fluxsign: '+' or '-' to specify which diretion normal to the surface for flux calculation
        sample_spacing: distance between sampled points for line plots
    meshparams = buildgeom.MeshParameters instance
  No return value.
  Output files are created."""

  #Mesh input files
  mesh_xml, surface_xml, volume_xml = solver_general.List_Mesh_Input_Files(modelparams.meshname)

  #Load mesh and meshfunctions
  mesh=Mesh(mesh_xml)
  surfaces=MeshFunction("size_t", mesh, surface_xml) #Mesh Function of Physical Surface number
  volumes=MeshFunction("size_t", mesh, volume_xml) #Mesh function of Physical Volume number

  #Function space for scalars and vectors
  V = FunctionSpace(mesh,'CG',1) #CG="continuous galerkin", ie "Lagrange"
  V_vec = VectorFunctionSpace(mesh, "CG", 1)

  #Dirichlet boundary conditions
  bcs=BCParameters(**modelparams.boundaryconditions).to_bclist()
  
  #Neumann boundary conditions
  #they are all zero in this case
  ds = Measure("ds",domain=mesh,subdomain_data=surfaces) ##TODO: specify which surfaces are Neumann?

  #Define variational problem
  c=TrialFunction(V)
  v=TestFunction(V)
  a=dot(grad(c),grad(v))*dx
  L=Constant(0)*v*ds
  
  #Solve
  c=Function(V)
  solve(a==L,c,bcs)

  #Output
  solver_general.Create_Output(modelparams,meshparams,c,modelparams.dataextraction)

  #Done
  return

#Support command-line arguments
if __name__ == '__main__':
  ##TODO: this needs to be updated
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Solve the unhomogenized fickian diffusion equation with fenics')
  parser.add_argument('bc_params_yaml', help='path to boundary conditions parameter yaml file')
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.bc_params_yaml), "Boundary conditions parameter definition file does not exist: %s"%cmdline.bc_params_yaml

  #Read in the yaml file
  solruns=useful.readyaml_multidoc(cmdline.bc_params_yaml)
  
  #Run each requested analysis
  for run in solruns:
    params=argparse.Namespace(**run)
    SolveMesh(params)

