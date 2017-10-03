#Solve the unhomogenized Fickian diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import os
import os.path as osp
import pickle
import sys

#Site packages
from fenics import *

#Local
sys.path.append(osp.abspath('..'))
import useful

#Constants
xmldir=osp.abspath('../mesh/xml')
solndir=osp.abspath('../solutions')
pickle_protocol = 4 #The newest protocol, requires python 3.4 or above.

def SolveMesh(basename,params):
  """Solve the unhomogenized Fickian diffusion equation on the indicated mesh.
  Arguments:
    basename = base name used for mesh, as string
    params = Namespace containing necessary parameter values:
      topsurf = physical surface number for top surface
      basesurf = physical surface number for base surface
      topval = concentration value at top surface
      baseval = concentration value at base surface
  No return value.
  Output files are created."""

  #Output location(s)
  outdir=osp.join(solndir,basename)
  if not osp.isdir(outdir):
    os.mkdir(outdir)
  pklfile=osp.join(outdir,'results.pkl')

  #Mesh input files
  mesh_xml=osp.join(xmldir,basename+'.xml')
  surface_xml=osp.join(xmldir,basename+'_facet_region.xml')
  volume_xml=osp.join(xmldir,basename+'_physical_region.xml')

  #Load mesh and meshfunctions
  mesh=Mesh(mesh_xml)
  surfaces=MeshFunction("size_t", mesh, surface_xml) #Mesh Function of Physical Surface number
  volumes=MeshFunction("size_t", mesh, volume_xml) #Mesh function of Physical Volume number

  #Function space for scalars and vectors
  V = FunctionSpace(mesh,'CG',1) #CG="continuous galerkin", ie "Lagrange"
  V_vec = VectorFunctionSpace(mesh, "CG", 1)

  #Dirichlet boundary conditions
  bcs=[]
  dpairs=[(params.basesurf,params.baseval), (params.topsurf,params.topval)] #Physical surface and Dirichlet value pairs
  for psurf,val in dpairs:
    bcs.append(DirichletBC(V,val,surfaces,psurf))
  
  #Neumann boundary conditions
  #they are all zero in this case
  ds = Measure("ds",domain=mesh,subdomain_data=boundaries) ##TODO: specify which surfaces are Neumann?

  #Define variational problem
  c=TrialFunction(V)
  v=TestFunction(V)
  a=dot(grad(c),grad(v))*dx
  L=Constant(0)*v*ds
  
  #Solve
  c=Function(V)
  solve(a==L,c,bcs)

  #VTK of concnentration
  vtk_c=File(osp.join(outdir,'conc.pvd'))
  vtk_c << c

  #Gradient
  j=project(Constant(-1)*grad(c),V_vec)
  vtk_j=File(osp.join(outdir,'flux.pvd'))
  vtk_j << j

  #Flux integral over surface
  ##TODO
  
  #Data for plots
  ##TODO

  #Pickle
  pobj={}
  ll=locals()
  for var in ['basename','params','mesh','surfaces','volumes','c','V','V_vec']:
    pobj[var]=ll[var]
  with open(pklfile,'wb') as fp:
    pickle.dump(pobj,fp,pickle_protocol)

  #Done
  return

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Solve the unhomogenized fickian diffusion equation with fenics')
  parser.add_argument('basename', help="base name of mesh file")
  parser.add_argument('paramfile', help="path to input parameters file")
  cmdline=parser.parse_args()

  #Read in the yaml file(s)
  paramsdict=useful.readyaml(cmdline.paramfile)
  params=argparse.Namespace(**paramsdict)

  #Create the file
  SolveMesh(cmdline.basename,params)
