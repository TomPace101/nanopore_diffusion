#Solve the unhomogenized Fickian diffusion problem
#and extract the data needed for post-processing efforts

#Standard library
import argparse
import os
import os.path as osp
import sys

#Site packages
from fenics import *
import numpy as np

#Local
import solver_general
import useful
import plotdata
from folderstructure import *

#TODO: there are probably parts of this that should be refactored into functions in a more general file, once one exists
def SolveMesh(params):
  """Solve the unhomogenized Fickian diffusion equation on the indicated mesh.
  Arguments:
    params = Namespace containing necessary parameter values:
      modelname = base name used for model, as string
      meshparams = name of mesh parameters file, as string
      meshname = base name used for mesh, as string
      topsurf = physical surface number for top surface
      basesurf = physical surface number for base surface
      topval = concentration value at top surface
      baseval = concentration value at base surface
      fluxsurf: physical surface number for flux measurement
      fluxsign: '+' or '-' to specify which diretion normal to the surface for flux calculation
      sample_spacing: distance between sampled points for line plots
  No return value.
  Output files are created."""

  #Read mesh parameters
  #TODO: it would be better if we didn't have to search through the whole file
  meshparamsfile=osp.join(srcfolder,params.meshparams) #TODO: it would be better if the paths was relative to params/mesh instead
  for meshobj in useful.readyaml_multidoc(meshparamsfile):
    if meshobj['meshname']==params.meshname:
      break
  meshparams=argparse.Namespace(**meshobj)

  #Output location(s)
  outdir=osp.join(solnfolder,params.modelname)
  if not osp.isdir(outdir):
    os.mkdir(outdir)
  pklfile=osp.join(outdir,'results.pkl')

  #Mesh input files
  mesh_xml, surface_xml, volume_xml = solver_general.List_Mesh_Input_Files(params)

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
  ds = Measure("ds",domain=mesh,subdomain_data=surfaces) ##TODO: specify which surfaces are Neumann?

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
  n=FacetNormal(mesh)
  dsi=Measure('interior_facet',domain=mesh,subdomain_data=surfaces)
  totflux=assemble(dot(j,n(params.fluxsign))*dsi(params.fluxsurf))

  #Effective Diffusion Constant
  cell_area = meshparams.Lx * meshparams.Ly
  c_samples=[c(0,0,zv) for zv in [meshparams.H, meshparams.H + meshparams.tm]]
  delta_c=c_samples[1]-c_samples[0]
  Deff=float(totflux/cell_area*meshparams.tm/delta_c)

  #Data for plots
  #Hard-coded centerline plot for now
  zr=np.arange(0, 2*meshparams.H + meshparams.tm, params.sample_spacing)
  zlist=[]
  clist=[]
  for z in zr:
    zlist.append(z)
    tup=(0,0,z)
    clist.append(c(*tup))
  zarr=np.array(zlist)
  carr=np.array(clist)
  pd=plotdata.PlotSeries(xvals=zarr,yvals=carr,label='concentration along centerline')
  pklfile=osp.join(outdir,'plotdata_CL_c.pkl')
  pd.to_pickle(pklfile)

  #Pickle
  #Nice try, but "can't pickle SwigPyOjbect objects"
  # pobj={}
  # ll=locals()
  # for var in ['params.meshname','params','mesh','surfaces','volumes','c','V','V_vec']:
  #   pobj[var]=ll[var]
  # useful.writepickle(pobj,pklfile)
  
  #Results yaml
  volfrac = np.pi*meshparams.R**2/(4*meshparams.Lx*meshparams.Ly)
  robj={'totflux':totflux, 'Deff':Deff, 'free_volume_frac':volfrac}
  robj.update(meshparams.__dict__)
  robj.update(params.__dict__)
  useful.writeyaml(robj,osp.join(outdir,'results.yaml'))

  #Done
  return

#Support command-line arguments
if __name__ == '__main__':
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

