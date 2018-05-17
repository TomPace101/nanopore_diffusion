
#Output for time-domain problems

#Standard library
import os.path as osp

#Site packages
import fenics as fem
import numpy as np

#Local
import folderstructure as FS
import plotdata


def calc_netcharge(self,attrname='netcharge',solnname='soln'):
  """Calculate the net charge distribution in the mesh
  Arguments:
    attrname = name of simulator attribute to store resulting field
    solnname = name of attribute storing the solution field as a list (all species and potential, probably created by splitfield)
  Required attributes:
  The new attribute is added.
  No return value.
  No output files created."""
  clist = getattr(self,solnname)
  rho_list=[]
  for s in range(self.Nspecies):
    rho_list.append(fem.Constant(self.species.z[s])*clist[s])
  rho = fem.project(sum(rho_list),mesh=self.meshinfo.mesh) #V=self.V.sub(0).collapse() would also work
  rho.rename('rho','charge density') #Without this, the function gets a different name each time it is created, which crashes paraview
  setattr(self,attrname,rho)
  return

def td_solutionfield(self,filename,attrname='soln',idx=None):
  """Write solution field to VTK file at each timestep
  Arguments:
    filename = name of output file, as string
      File will be created in the output directory (self.outdir)
    attrname = name of attribute to output, as string, defaults to 'soln'
    idx = index of the solution field to write out, None (default) if not a sequence
  Required attributes:
    outdir = output directory, as string
    soln (or other given by attrname) = FEniCS Function containing solution
  New attributes:
    td_vtk_files = dictionary of FEniCS File objects, by filename
  No return value.
  Output file is created on first call, and updated each subsequent call."""
  #Create dictionary of files if not yet present
  if not hasattr(self,'td_vtk_files'):
    self.td_vtk_files={}
  #Create new file if not already present
  if not filename in self.td_vtk_files.keys():
    self.td_vtk_files[filename]=fem.File(osp.join(self.outdir,filename))
  #Output data
  output = getattr(self,attrname)
  if idx is not None:
    output = output[idx]
  self.td_vtk_files[filename] << (output,self.t)
  return

def td_pointhistory(self,location,plotname,label,attrname='soln',idx=None):
  """Get solution value at a single model point at each timestep
  Arguments:
    specifier = argument to get_pointcoords
    plotname = name of plot in outdata.plots, as string
    label = series label to assign, as string
    attrname = name of attribute to output, as string, defaults to 'soln'
  Required attributes:
    outdata = instance of OutData
    parametric_locations = only required if needed by location specifier, dictionary of parametric locations
  New attributes:
    td_point_series = dictionary of plotdata.Series objects
  No return value.
  No output file generated."""
  #Create dictionary of point histories if not yet present
  if not hasattr(self,'td_point_series'):
    self.td_point_series={}
  #Create new series if not already present
  if not (plotname,label) in self.td_point_series.keys():
    newseries=plotdata.PlotSeries(xvals=np.array([]),yvals=np.array([]),label=label,metadata={})
    self.td_point_series[(plotname,label)]=newseries
    if plotname not in self.outdata.plots.keys():
      self.outdata.plots[plotname]=[]
    self.outdata.plots[plotname].append(newseries)
  #Get coordinates of specified point
  coords=self.get_pointcoords(location)
  #Store data
  datafunction=getattr(self,attrname)
  if idx is not None:
    datafunction = datafunction[idx]
  outval=datafunction(*coords)
  series=self.td_point_series[(plotname,label)]
  series.xvals=np.append(series.xvals,self.t)
  series.yvals=np.append(series.yvals,outval)

