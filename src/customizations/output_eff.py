
#Output functions originally for 3D pore problems

#Standard library
import os.path as osp

#Site packages
import fenics as fem
import numpy as np

#Local
import folderstructure as FS
import plotdata

def fluxfield(self,filename):
  """Flux as vector field (new attribute, and VTK file)
  Arguments:
    filename = name of output file, as string
      File will be created in the output directory (self.outdir)
  Required attributes:
    conditions.D_bulk = bulk diffusion constant for the medium
    outdir = output directory, as string
    soln = FEniCS Function containing solution
    V_vec = FEniCS VectorFunctionSpace
  New attributes:
    flux = flux field as a MeshFunction
  No return value.
  Output file is written."""
  D_bulk=self.conditions.D_bulk
  self.flux=fem.project(fem.Constant(-D_bulk)*fem.grad(self.soln),self.V_vec)
  vtk_file=fem.File(osp.join(self.outdir,filename))
  vtk_file << self.flux
  return

def fluxintegral(self,fluxsurf,name,internal=False,fluxsign=None,normalvar=None): ##TODO: store also quadrupled value for unit cell?
  """Flux integral over specified facet
  Arguments:
    fluxsurf = physical facet number for flux measurement
    name = name for storage in the results dictionary
    internal = boolean, default False, True to use internal boundary, False for external
    fluxsign = '+' or '-' to specify which direction normal to the facet for flux calculation
      Required only if internal==True
    normalvar = optional variable name to write the facet normal components to, as a sequence
  Required attributes:
    flux = flux as vector field
      This requires a previous call to fluxfield
    mesh = FEniCS Mesh object
    facet = FEniCS MeshFunction object for facet numbers
  No new attributes.
  New item(s) added to results dictionary.
  No return value.
  No output files."""
  n=fem.FacetNormal(self.mesh)
  if internal:
    integral_type='interior_facet'
    assert fluxsign=='+' or fluxsign=='-', "Invalid fluxsign: %s"%str(fluxsign)
    this_n=n(fluxsign)
  else:
    integral_type='exterior_facet'
    this_n=n
  if normalvar is not None:
    self.results[normalvar]=['not_yet_computed'] ##TODO: find a way to get coordinates of the facet normal
  this_ds=fem.Measure(integral_type,domain=self.mesh,subdomain_data=self.facets)
  totflux=fem.assemble(fem.dot(self.flux,this_n)*this_ds(fluxsurf))
  self.results[name]=totflux
  return

def effective_diffusion(self,name,totflux_name):
  """Calculate effective diffusion constant
  Arguments:
    name = name for storage in the results dictionary
    totflux_name = name of previously calculated total flux in results dictionary
      This requires a previous call to fluxintegral.
  Required attributes:
    tmplvalues = tmplvalues attribute of buildgeom.MeshParameters object
    results[toflux_name] = result from previous call to fluxintegral
  No new attributes.
  New item added to results dictionary.
  No return value.
  No output files."""
  quarter_area = self.tmplvalues['Lx'] * self.tmplvalues['Ly']
  samples=[self.soln(0,0,zv) for zv in [self.tmplvalues['H'], self.tmplvalues['H'] + self.tmplvalues['tm']]]
  delta=samples[1]-samples[0]
  Deff=float(self.results[totflux_name]/quarter_area*self.tmplvalues['tm']/delta)
  self.results[name]=Deff
  return

def volfrac(self,name):
  """Calculate free volume fraction
  Arguments:
    name = name for storage in the results dictionary
  Required attributes:
    tmplvalues = tmplvalues attribute of buildgeom.MeshParameters object
  No new attributes.
  New item added to results dictionary.
  No return value.
  No output files."""
  self.results[name]=np.pi*self.tmplvalues['R']**2/(4*self.tmplvalues['Lx']*self.tmplvalues['Ly'])
  return

#We could generalize this by specifying:
# - a center point
# - an axis of rotation defining the plane
# - a point to define the orientation of theta=0
def profile_radial(self,zval,theta,num,plotname,label,attrname='soln',idx=None):
  """Data for plot of solution along radial line at specified Z, in specified direction
  Arguments:
    zval = 
    theta = theta-angle in degrees from x-axis, as float
    num = number of sampled points
    indep = index of the coordinate parameter to use as the independent variable for the plot (zero-based)
    plotname = name of plot in outdata.plots, as string
    label = series label to assign, as string
    attrname = name of attribute to output, as string, defaults to 'soln'
  Required attributes:
    outdata = instance of OutData
    mesh_metadata = dictionary of mesh metadata
      MUST contain radius under key 'R'
  No new attributes.
  Nothing added to results dictionary.
  No return value.
  Series is added to `outdata.plots`."""
  #Calculate start and end locations
  startloc=(0,0,zval)
  rads=np.radians(theta)
  xend=self.mesh_metadata['R']*np.cos(rads)
  yend=self.mesh_metadata['R']*np.sin(rads)
  endloc=(xend,yend,zval)
  #Call line_profile to finish
  self.line_profile(startloc,endloc,num,plotname,label,attrname,None,idx)

def OLD_profile_radial(self,spacing,plotname,label,theta,attrname='soln'):
  """Data for plot of solution along radial line at model mid-height, in specified direction
  Arguments:
    spacing = distance between sampled points for line plots
    plotname = name of plot in outdata.plots, as string
    label = series label to assign, as string
    theta = theta-angle in degrees from x-axis, as float
    attrname = name of attribute to output, as string, defaults to 'soln'
  Required attributes:
    outdata = instance of OutData
    tmplvalues = tmplvalues attribute of buildgeom.MeshParameters object
    modelparams = solver_run.ModelParameters object
    outdir = path to output directory, as string
  No new attributes.
  Nothing added to results dictionary.
  No return value.
  Series is added to `outdata.plots`."""
  #Get the object with the data
  vals=getattr(self,attrname)
  #Extract data points
  zval=self.tmplvalues['H'] + self.tmplvalues['tm']/2 #mid-height
  rads=np.radians(theta)
  cos=np.cos(rads)
  sin=np.sin(rads)
  tree=self.mesh.bounding_box_tree()
  rrange=np.arange(0,self.tmplvalues['R']+spacing,spacing)
  rlist=[]
  vlist=[]
  tuplist=[]
  for r in rrange:
    xval=r*cos
    yval=r*sin
    # if direction=='x':
    #   xval=r
    #   yval=0
    # else:
    #   xval=0
    #   yval=r
    tup=(xval,yval,zval)
    pt=fem.Point(*tup)
    inside=tree.collides(pt)
    #Is this point inside the mesh (including on the boundary)
    if inside:
      rlist.append(r)
      try:
        vlist.append(vals(*tup))
        inside='FalsePositive'
      except RuntimeError:
        #Couldn't get value for this point, so drop it
        rlist.pop(-1)
    #Track result
    tuplist.append((tup,inside))
  #Create PlotSeries
  rarr=np.array(rlist)
  varr=np.array(vlist)
  meta={'tuplist':tuplist}
  series=plotdata.PlotSeries(xvals=rarr,yvals=varr,label=label,metadata=meta)
  #Store data
  if not plotname in self.outdata.plots.keys():
    self.outdata.plots[plotname]=[]
  self.outdata.plots[plotname].append(series)
  return
