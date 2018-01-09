#Functions used by various solvers

#Standard library
import os
import os.path as osp
import sys

#Site packages
import fenics as fem
import numpy as np

#Local
from folderstructure import *
import useful
import buildgeom
import plotdata

#fenics version check
target_fenics_version='2017.1.0'
fenics_version_msg_template="This code was written for FEniCS version '%s'. Detected version '%s'. You can try to run it by uncommenting this line, but it may not work, due to FEniCS API changes."
assert fem.DOLFIN_VERSION_STRING == target_fenics_version, fenics_version_msg_template%(target_fenics_version,fem.DOLFIN_VERSION_STRING)

class GenericConditions(useful.ParameterSet):
  """Condition defnitions, to be subclassed by each equation as needed
  Attributes:
    elementorder = integer specifying equation order (1=1st, 2=2nd, etc) for finite elements
    bcdict = dictionary of Dirichlet boundary conditions: {physical surface number: solution value, ...}"""
  __slots__=['elementorder','bcdict']

class OutData(useful.ParameterSet):
  """Data from solver to be written to disk.
  Attributes:
    plots = Dictionary of plot data, {plotname: [plotdata.PlotSeries, ...], ...}"""
  __slots__=['plots']

class ModelParameters(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store generic solver parameters
  Attributes:
    modelname = stem name for output files
    meshparamsfile = name of yaml file containing the mesh defined by meshname (include .yaml extension)
    meshname = stem name for mesh files
    equation = name of equation to be solved
    properties = dictionary of property values
    conditions = parameters specifying boundary conditions, initial conditions, property values, etc.
      The parameters specified are specific to the equation being solved.
      See the equation's module for the equation-specific conditions object,
      which will usually inherit from GenericConditions.
    dataextraction = a sequence of data extraction commands
      Each command is a pair (cmdname, arguments), where
        cmdname = name of the solver object's method to call, as a string
        arguments = dictionary of arguments to the method: {argname: value,...}"""
  __slots__=('modelname','meshparamsfile','meshname','equation','conditions','dataextraction')

def List_Mesh_Input_Files(meshname,basedir):
  mesh_xml=osp.join(xmlfolder,basedir,meshname+'.xml')
  surface_xml=osp.join(xmlfolder,basedir,meshname+'_facet_region.xml')
  volume_xml=osp.join(xmlfolder,basedir,meshname+'_physical_region.xml')
  return mesh_xml, surface_xml, volume_xml

loaded_meshfiles=[]
loaded_meshes={}
def complete_by_ModelParameters(modelparams,solverclasses):
  """Run a model, given only its modelparameters object.
  Arguments:
    modelparams=ModelParameters instance
    solverclasses = dictionary of solver classes: {equation name: instance of a subclass of GenericSolver}"""
  #Load the mesh file if not already loaded
  if not modelparams.meshparamsfile in loaded_meshfiles:
    meshparams_fpath=osp.join(params_mesh_folder,modelparams.meshparamsfile)
    meshparams_gen=buildgeom.MeshParameters.all_from_yaml(meshparams_fpath)
    for mp in meshparams_gen:
      loaded_meshes[mp.meshname]=mp
    loaded_meshfiles.append(modelparams.meshparamsfile)
  #Get the requested mesh
  meshparams=loaded_meshes[modelparams.meshname]
  #Run the solver
  solverclasses[modelparams.equation].complete(modelparams,meshparams)

class GenericSolver:
  """A generic solver, to be subclassed by solvers for the specific equations
  This class is not directly usable itself.
  Derived classes should, at a minimum:
    - override __init__ to set up the variational problem (it's ok to use super() to set up the mesh)
    - implement a "solve" method to generate the solution
  and other data needed by their data extraction functions.
  Subclasses may choose to override the extraction functions provided here.
  Attributes:
    modelparams = ModelParameters instance
    meshparams = buildgeom.MeshParameters instance
    mesh = FEniCS Mesh
    surfaces = FEniCS MeshFunction of Physical Surface nubmer
    volumes = FEniCS MeshFunction of Physical Volume number
  """
  def __init__(self,modelparams,meshparams):
    """Initialize the solver by loading the Mesh and MeshFunctions.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
    #Store defining ParameterSet objects
    self.modelparams=modelparams
    self.meshparams=meshparams

    #Initialize output attributes
    self.results={}
    self.info=self.modelparams.to_dict()
    self.info['meshparams']=self.meshparams.to_dict()
    self.outdata=OutData(plots={})

  def loadmesh(self):
    """Load the mesh from the usual file locations"""
    #Mesh input files
    mesh_xml, surface_xml, volume_xml = List_Mesh_Input_Files(self.modelparams.meshname,self.meshparams.basename)

    #Load mesh and meshfunctions
    self.mesh=fem.Mesh(mesh_xml)
    self.surfaces=fem.MeshFunction("size_t", self.mesh, surface_xml) #Mesh Function of Physical Surface number
    self.volumes=fem.MeshFunction("size_t", self.mesh, volume_xml) #Mesh function of Physical Volume number

    return

  @classmethod
  def complete(cls,*args,diskwrite=True,as_action=False):
    """Convenience function to set up and solve the model, then generate all the requested output.
    Arguments:
      *args to be passed to the the solver class __init__
      diskwrite = boolean, True to write info and output data to disk
      as_action = boolean, True to return None rather than the solver
        This is used to allow this method to be used as a task in doit,
        which restricts what a task function can return.
    In general, it would be a bad idea to use diskwrite=False and as_action=True,
    because you'd then have no way to get any results from the solution at all."""
    obj=cls(*args)
    obj.solve()
    obj.create_output(diskwrite)
    if as_action:
      obj = None
    return obj

  def solve(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'solve' method."%str(type(self)))

  def create_output(self,diskwrite=True):
    """Process the data extraction commands
    Arguments:
      diskwrite = boolean, optional, True to write the output files.
        If false, the output will only be stored in the object itself
    Adds the following attributes:
      outdir = path to output directory, as string
      results = dictionary of input and output values
    Output files are generated."""
    #Output location(s)
    self.outdir=osp.join(solnfolder,self.modelparams.basename,self.modelparams.modelname)
    if not osp.isdir(self.outdir):
      os.makedirs(self.outdir)

    #Process each command
    for cmd in getattr(self.modelparams,'dataextraction',[]):
      #Function name and arguments
      funcname, kwargs = cmd
      #Call it
      try:
        getattr(self,funcname)(**kwargs)
      except Exception as einst:
        print("Excption occured for command: %s"%str(cmd), file=sys.stderr)
        raise einst

    #Put results into info
    self.info['results']=self.results
    #Write output files if requested
    if diskwrite:
      useful.writeyaml(self.info,osp.join(self.outdir,'info.yaml'))
      self.outdata.to_pickle(osp.join(self.outdir,'outdata.pkl'))

    #Done
    return

  def solutionfield(self,filename,attrname='soln'):
    """Write solution field to VTK file
    Arguments:
      filename = name of output file, as string
        File will be created in the output directory (self.outdir)
      attrname = name of attribute to output, as string, defaults to 'soln'
    Required attributes:
      outdir = output directory, as string
      soln = FEniCS Function containing solution
    No new attributes.
    No return value.
    Output file is written."""
    vtk_file=fem.File(osp.join(self.outdir,filename))
    vtk_file << getattr(self,attrname)
    return

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
    """Flux integral over specified surface
    Arguments:
      fluxsurf = physical surface number for flux measurement
      name = name for storage in the results dictionary
      internal = boolean, default False, True to use internal boundary, False for external
      fluxsign = '+' or '-' to specify which diretion normal to the surface for flux calculation
        Required only if internal==True
      normalvar = optional variable name to write the surface normal components to, as a sequence
    Required attributes:
      flux = flux as vector field
        This requires a previous call to fluxfield
      mesh = FEniCS Mesh object
      surface = FEniCS MeshFunction object for surface numbers
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
      self.results[normalvar]=['not_yet_computed'] ##TODO: find a way to get coordinates of the surface normal
    this_ds=fem.Measure(integral_type,domain=self.mesh,subdomain_data=self.surfaces)
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
      meshparams = buildgeom.MeshParameters object
      results[toflux_name] = result from previous call to fluxintegral
    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    quarter_area = self.meshparams.Lx * self.meshparams.Ly
    samples=[self.soln(0,0,zv) for zv in [self.meshparams.H, self.meshparams.H + self.meshparams.tm]]
    delta=samples[1]-samples[0]
    Deff=float(self.results[totflux_name]/quarter_area*self.meshparams.tm/delta)
    self.results[name]=Deff
    return

  def volfrac(self,name):
    """Calculate free volume fraction
    Arguments:
      name = name for storage in the results dictionary
    Required attributes:
      meshparams = buildgeom.MeshParameters object
    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    self.results[name]=np.pi*self.meshparams.R**2/(4*self.meshparams.Lx*self.meshparams.Ly)
    return

  def profile_centerline(self,spacing,plotname,label,attrname='soln'):
    """Data for plot of solution profile along centerline
    Arguments:
      spacing = distance between sampled points for line plots
      plotname = name of plot in outdata.plots, as string
      label = series label to assign, as string
      attrname = name of attribute to output, as string, defaults to 'soln'
    Required attributes:
      meshparams = buildgeom.MeshParameters object
      modelparams = ModelParameters object
      outdir = path to output directory, as string
    No new attributes.
    Nothing added to results dictionary.
    No return value.
    Series is added to `outdata.plots`."""
    #Get the object with the data
    vals=getattr(self,attrname)
    #Extract data points
    zr=np.arange(0, 2*self.meshparams.H + self.meshparams.tm, spacing)
    zlist=[]
    vlist=[]
    for z in zr:
      zlist.append(z)
      tup=(0,0,z)
      vlist.append(vals(*tup))
    #Create PlotSeries
    zarr=np.array(zlist)
    varr=np.array(vlist)
    series=plotdata.PlotSeries(xvals=zarr,yvals=varr,label=label)
    #Store data
    if not plotname in self.outdata.plots.keys():
      self.outdata.plots[plotname]=[]
    self.outdata.plots[plotname].append(series)
    return

  def profile_radial(self,spacing,plotname,label,theta,attrname='soln'):
    """Data for plot of solution along radial line at model mid-height, in specified direction
    Arguments:
      spacing = distance between sampled points for line plots
      plotname = name of plot in outdata.plots, as string
      label = series label to assign, as string
      theta = theta-angle in degrees from x-axis, as float
      attrname = name of attribute to output, as string, defaults to 'soln'
    Required attributes:
      meshparams = buildgeom.MeshParameters object
      modelparams = ModelParameters object
      outdir = path to output directory, as string
    No new attributes.
    Nothing added to results dictionary.
    No return value.
    Series is added to `outdata.plots`."""
    #Get the object with the data
    vals=getattr(self,attrname)
    #Extract data points
    zval=self.meshparams.H + self.meshparams.tm/2 #mid-height
    rads=np.radians(theta)
    cos=np.cos(rads)
    sin=np.sin(rads)
    tree=self.mesh.bounding_box_tree()
    rrange=np.arange(0,self.meshparams.R+spacing,spacing)
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
