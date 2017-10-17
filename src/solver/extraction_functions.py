#Functions called by solver_general.Create_Output

#Standard library
import sys
import os.path as osp

#Site packages
from fenics import *
import numpy as np

#Local
if not osp.abspath('..') in sys.path:
  sys.path.append(osp.abspath('..'))
from folderstructure import *
import plotdata


#All functions must accept the following arguments:
# solution,results,outdir,**kwargs
#No return values are used.

def solutionfield(soln,results,outdir,filename):
  "Write solution field to VTK file"
  vtk_file=File(osp.join(outdir,filename))
  vtk_file << soln
  return

def fluxfield(soln,results,outdir,filename):
  "Flux as vector field"
  flux=project(Constant(-1)*grad(soln),V_vec) ##TODO: don't have V_vec
  vtk_file=File(osp.join(outdir,filename))
  vtk_file << flux
  return

def fluxintegral(soln,results,outdir,fluxsurf,fluxsign,name):
  "Flux integral over specified surface"
  n=FacetNormal(mesh) ##TODO: don't have mesh
  dsi=Measure('interior_facet',domain=mesh,subdomain_data=surfaces) ##TODO: don't have mesh or surfaces
  totflux=assemble(dot(j,n(fluxsign))*dsi(fluxsurf))
  results[name]=totflux
  return

def effective_diffusion(soln,results,outdir,name):
  "Calculate effective diffusion constant"
  quarter_area = meshparams.Lx * meshparams.Ly ##TODO: don't have meshparams
  samples=[soln(0,0,zv) for zv in [meshparams.H, meshparams.H + meshparams.tm]] ##TODO: don't have meshparams
  delta=samples[1]-samples[0]
  Deff=float(totflux/quarter_area*meshparams.tm/delta) #TODO: don't have meshparams, and need total flux (order of operations matters!)
  results[name]=Deff
  return

def volfrac(soln,results,outdir,name):
  "Calculate free volume fraction"
  volfrac = np.pi*meshparams.R**2/(4*meshparams.Lx*meshparams.Ly) ##TODO: don't have meshparams
  results[name]=volfrac
  return

def profile_centerline(soln,results,outdir,spacing,filename,label):
  "Data for plot of concentration profile along centerline"
  ##TODO: replace c with soln
  zr=np.arange(0, 2*meshparams.H + meshparams.tm, spacing) ##TODO: don't have meshparams
  zlist=[]
  vlist=[]
  for z in zr:
    zlist.append(z)
    tup=(0,0,z)
    vlist.append(soln(*tup))
  zarr=np.array(zlist)
  varr=np.array(vlist)
  meta=dict([(k,getattr(meshparams,k)) for k in ['Lx','Ly','R','tm','H']]) ##TODO: don't have meshparams
  meta.update(bcvals.__dict__) #TODO: don't have bcvals
  pd=plotdata.PlotSeries(xvals=zarr,yvals=varr,label=label,metadata=meta)
  pklfile=osp.join(outdir,filename)
  pd.to_pickle(pklfile)
  return

fdir=globals()
exfunc_names=['solutionfield','fluxfield','fluxintegral','effective_diffusion','volfrac','profile_centerline']
exfuncs=dict([(fn,fdir[fn]) for fn in exfunc_names])

