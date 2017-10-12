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
from folderstructure import *
import solver_general
import useful
import plotdata

#TODO: there are probably parts of this that should be refactored into functions in a more general file, once one exists
def SolveMesh(modelparams, meshparams):
  """Solve the unhomogenized Fickian diffusion equation on the indicated mesh.
  Arguments:
    modelparams = solver_general.ModelParameters instance
      boundaryconditions:
        topsurf = physical surface number for top surface
        basesurf = physical surface number for base surface
        topval = concentration value at top surface
        baseval = concentration value at base surface
      dataextraction:
        fluxsurf: physical surface number for flux measurement
        fluxsign: '+' or '-' to specify which diretion normal to the surface for flux calculation
        sample_spacing: distance between sampled points for line plots
    meshparams = buildgeom.MeshParameters instance
  No return value.
  Output files are created."""

  #Output location(s)
  outdir=osp.join(solnfolder,modelparams.modelname)
  if not osp.isdir(outdir):
    os.mkdir(outdir)

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
  bcs=[]
  bcvals=argparse.Namespace(**modelparams.boundaryconditions)
  dpairs=[(bcvals.basesurf,bcvals.baseval), (bcvals.topsurf,bcvals.topval)] #Physical surface and Dirichlet value pairs
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

  #Data extraction
  extraction=argparse.Namespace(**modelparams.dataextraction)

  #Flux integral over surface
  n=FacetNormal(mesh)
  dsi=Measure('interior_facet',domain=mesh,subdomain_data=surfaces)
  totflux=assemble(dot(j,n(extraction.fluxsign))*dsi(extraction.fluxsurf))

  #Effective Diffusion Constant
  cell_area = meshparams.Lx * meshparams.Ly
  c_samples=[c(0,0,zv) for zv in [meshparams.H, meshparams.H + meshparams.tm]]
  delta_c=c_samples[1]-c_samples[0]
  Deff=float(totflux/cell_area*meshparams.tm/delta_c)

  #Data for plots
  #Hard-coded centerline plot for now
  zr=np.arange(0, 2*meshparams.H + meshparams.tm, extraction.sample_spacing)
  zlist=[]
  clist=[]
  for z in zr:
    zlist.append(z)
    tup=(0,0,z)
    clist.append(c(*tup))
  zarr=np.array(zlist)
  carr=np.array(clist)
  meta=dict([(k,getattr(meshparams,k)) for k in ['Lx','Ly','R','tm','H']])
  meta.update(bcvals.__dict__)
  pd=plotdata.PlotSeries(xvals=zarr,yvals=carr,label='concentration along centerline',metadata=meta)
  pklfile=osp.join(outdir,'plotdata_CL_c.pkl')
  pd.to_pickle(pklfile)

  #Pickle
  #Nice try, but "can't pickle SwigPyOjbect objects"
  # pobj={}
  # ll=locals()
  # for var in ['modelparams','mesh','surfaces','volumes','c','V','V_vec']:
  #   pobj[var]=ll[var]
  # useful.writepickle(pobj,pklfile)
  
  #Results yaml
  volfrac = np.pi*meshparams.R**2/(4*meshparams.Lx*meshparams.Ly)
  robj={'totflux':totflux, 'Deff':Deff, 'free_volume_frac':volfrac}
  robj.update(meshparams.to_dict())
  robj.update(modelparams.to_dict())
  useful.writeyaml(robj,osp.join(outdir,'results.yaml'))

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

