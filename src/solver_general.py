#Functions used by various solvers

#Standard library
import os
import os.path as osp
import sys

#Site packages
import fenics as fem
import numpy as np

#Local
import folderstructure as FS
import common
import buildgeom
import plotdata

#fenics version check
target_fenics_version='2017.1.0'
fenics_version_msg_template="This code was written for FEniCS version '%s'. Detected version '%s'. You can try to run it by uncommenting this line, but it may not work, due to FEniCS API changes."
assert fem.DOLFIN_VERSION_STRING == target_fenics_version, fenics_version_msg_template%(target_fenics_version,fem.DOLFIN_VERSION_STRING)

class GenericConditions(common.ParameterSet):
  """Condition defnitions, to be subclassed by each equation as needed
  Attributes:
    elementorder = integer specifying equation order (1=1st, 2=2nd, etc) for finite elements
    bcdict = dictionary of Dirichlet boundary conditions: {physical facet number: solution value, ...}"""
  __slots__=['elementorder','bcdict']

class OutData(common.ParameterSet):
  """Data from solver to be written to disk.
  Attributes:
    plots = Dictionary of plot data, {plotname: [plotdata.PlotSeries, ...], ...}"""
  __slots__=['plots']

class ModelParametersBase(common.ParameterSet):
  """Subclass of common.ParameterSet to store generic solver parameters
  Attributes:
    To be read in from file:
      modelname = stem name for output files
      meshparamsfile = name of yaml file containing the mesh defined by meshname (include .yaml extension)
      meshname = stem name for mesh files
      equation = name of equation to be solved
      properties = dictionary of property values
      conditions = parameters specifying boundary conditions, initial conditions, property values, etc.
        The parameters specified are specific to the equation being solved.
        See the equation's module for the equation-specific conditions object,
        which will usually inherit from GenericConditions.
      dataextraction = a sequence of data extraction commands to execute after solving the model
        Each command is a pair (cmdname, arguments), where
          cmdname = name of the solver object's method to call, as a string
          arguments = dictionary of arguments to the method: {argname: value,...}
      datasteps = a sequence of data extraction commands to execute at each time step
        Note: Not all solvers will support this argument. It is mainly intended for time-domain equations.
        The command structure is the same as for dataextraction.
    To be calculated from methods:
      outdir = folder containing output files"""
  __slots__=('modelname','meshparamsfile','meshname','equation','conditions','dataextraction','datasteps','outdir')
  _required_attrs=['modelname','meshparamsfile','meshname','equation','conditions']
  _config_attrs=['basename']+_required_attrs+['dataextraction','datasteps']
  _taskname_src_attr='modelname'

  def __init__(self,**kwd):
    #Initialization from base class
    super().__init__(**kwd)
    self.outdir=osp.join(FS.solnfolder,self.basename,self.modelname)

class GenericSolver:
  """A generic solver, to be subclassed by solvers for the specific equations
  This class is not directly usable itself.
  Derived classes should, at a minimum:
    - override __init__ to set up the variational problem (it's ok to use super() to set up the mesh)
    - implement a "solve" method to generate the solution
  and other data needed by their data extraction functions.
  Subclasses may choose to override the extraction functions provided here.
  Attributes:
    modelparams = solver_run.ModelParameters instance
    meshparams = buildgeom.MeshParameters instance
    tmplvalues = shortcut attribute for meshparams.tmplvalues
    mesh = FEniCS Mesh
    facets = FEniCS MeshFunction of gmsh Physical Surface number (3D) or Physical Line number (2D)
    cells = FEniCS MeshFunction of gmsh Physical Volume number (3D) or Physical Surface number (2D)
  """
  def __init__(self,modelparams,meshparams):
    """Initialize the solver by loading the Mesh and MeshFunctions.
    Arguments:
      modelparams = solver_run.ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
    #Store defining ParameterSet objects
    self.modelparams=modelparams
    self.meshparams=meshparams

    #Initialize output attributes and intermediates
    self.results={}
    self.info=self.modelparams.config_dict
    self.info['meshparams']=self.meshparams.config_dict
    self.outdata=OutData(plots={})
    self.tmplvalues=self.meshparams.tmplvalues

  def loadmesh(self):
    """Load the mesh from the usual file locations
    A note on the terminology used in FEniCS and gmsh:
      The FEniCS parts are from page 185-186 of the FEniCS book
      d = number of dimensions in entity, D = number of dimensions in problem (maximum entity dimension)
      D-d = "codimension" of entity
      Terms:
        D=2, d=1: fenics facet (facet_region xml) = fenics edge = gmsh physical line
        D=2, d=2: fenics cell (physical_region xml) = fenics face = gmsh physical surface
        D=3, d=2: fenics facet (facet_region xml) = fenics face = gmsh physical surface
        D=3, d=3: fenics cell (physical_region xml) = fenics ____ = gmsh physical volume
        also, d=0 is a fenics vertex
      """
    #Load mesh and meshfunctions
    self.mesh=fem.Mesh(self.modelparams.mesh_xml)
    self.facets=fem.MeshFunction("size_t", self.mesh, self.modelparams.facet_xml)
    self.cells=fem.MeshFunction("size_t", self.mesh, self.modelparams.cell_xml)
    return

  @classmethod
  def complete(cls,*args,diskwrite=True):
    """Convenience function to set up and solve the model, then generate all the requested output.
    Arguments:
      *args to be passed to the the solver class __init__
      diskwrite = boolean, True to write info and output data to disk
    In general, it would be a bad idea to use diskwrite=False and as_action=True,
    because you'd then have no way to get any results from the solution at all."""
    obj=cls(*args)
    obj.solve()
    obj.create_output(diskwrite)
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
    self.outdir=self.modelparams.outdir
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
      common.writeyaml(self.info,osp.join(self.outdir,'info.yaml'))
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
    """Flux integral over specified facet
    Arguments:
      fluxsurf = physical facet number for flux measurement
      name = name for storage in the results dictionary
      internal = boolean, default False, True to use internal boundary, False for external
      fluxsign = '+' or '-' to specify which diretion normal to the facet for flux calculation
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

  def profile_centerline(self,spacing,plotname,label,attrname='soln'):
    """Data for plot of solution profile along centerline
    Arguments:
      spacing = distance between sampled points for line plots
      plotname = name of plot in outdata.plots, as string
      label = series label to assign, as string
      attrname = name of attribute to output, as string, defaults to 'soln'
    Required attributes:
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
    zr=np.arange(0, 2*self.tmplvalues['H'] + self.tmplvalues['tm'], spacing)
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
