#Functions used by various solvers

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import os
import os.path as osp
import sys
import types

#Site packages
import fenics as fem
import numpy as np

#Local
import folderstructure as FS
import common
import plotdata

class GenericConditions(common.ParameterSet):
  """Condition defnitions, to be subclassed by each equation as needed
  Attributes:
    elementorder = integer specifying equation order (1=1st, 2=2nd, etc) for finite elements"""
  __slots__=['elementorder']

class OutData(common.ParameterSet):
  """Data from solver to be written to disk.
  Attributes:
    plots = Dictionary of plot data, {plotname: [plotdata.PlotSeries, ...], ...}"""
  __slots__=['plots']

class SolverCustomizations(common.ParameterSet):
  """Information about custom methods and attributes for the solver class
  Attributes:
    modules = a sequence of modules to be imported. All functions defined inside become methods of the solver.
      (Technically, all functions whose names appear in dir(module), which could be tailored by defining __dir__ if desired.)
    initializations = a dictionary {module name: {variable: value}}
      Upon loading the module, the values in this dictionary will be passed as keyword arguments
        to the function `initialize_module`, if present, within the module.
      Modules listed here but not in `modules` are silently ignored.
    extra = dictionary {additional solver attributes: initial value}"""
  __slots__=['modules','initializations','extra']

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
      customizations = parameters specifying an instance of SolverCustomizations
    To be calculated from methods:
      outdir = folder containing output files"""
  __slots__=('modelname','meshparamsfile','meshname','equation','conditions','dataextraction','datasteps','customizations','outdir')
  _required_attrs=['modelname','meshparamsfile','meshname','equation','conditions']
  _config_attrs=['basename']+_required_attrs+['dataextraction','datasteps']
  _taskname_src_attr='modelname'

  def __init__(self,**kwd):
    #Initialization from base class
    super(ModelParametersBase, self).__init__(**kwd)
    #Get output directory
    self.outdir=osp.join(FS.solnfolder,self.basename,self.modelname)
    #Load customization info, if present
    if hasattr(self,'customizations'):
      self.customizations=SolverCustomizations(**self.customizations)

class GenericSolver(object):
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
    self.diskwrite=True
    self.outdir=self.modelparams.outdir
    self.results={}
    self.info=self.modelparams.config_dict
    self.info['meshparams']=self.meshparams.config_dict
    self.outdata=OutData(plots={})
    self.tmplvalues=self.meshparams.tmplvalues
    
    #Apply customizations
    if hasattr(self.modelparams,'customizations'):
      loaded_modules={}
      module_initializations=getattr(self.modelparams.customizations,'initializations',{})
      #Bind methods
      for modname in getattr(self.modelparams.customizations,'modules',[]):
        #Load module
        themod=importlib.import_module(modname)
        loaded_modules[modname]=themod
        #Intialize, if requested
        if modname in module_initializations.keys():
          kwargs = module_initializations[modname]
          if (kwargs is not None) and hasattr(themod,'initialize_module'):
            themod.initialize_module(**kwargs)
        #Assign all module functions as solver methods
        mod_contents=dict([(f,getattr(themod,f)) for f in dir(themod)])
        for nm, itm in mod_contents.items():
          if isinstance(itm,types.FunctionType):
            setattr(self,nm,types.MethodType(itm,self)) #Must type cast to MethodType in order to get implicit first argument `self`
      #Assign extra attributes
      for k,v in getattr(self.modelparams.customizations,'extra',{}).items():
        setattr(self,k,v)

  def loadmesh(self):
    """Load the mesh from the usual file locations
    Also load the parametric locations file, if present.
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
    #Load parametric locations file, if found
    if osp.isfile(self.modelparams.paramlocsfile):
      self.parametric_locations=common.readyaml(self.modelparams.paramlocsfile)
    return

  @classmethod
  def complete(cls,*args):
    """Convenience function to set up and solve the model, then generate all the requested output.
    Arguments:
      *args to be passed to the the solver class __init__"""
    obj=cls(*args)
    obj.solve()
    obj.create_output()
    return obj

  def solve(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'solve' method."%str(type(self)))

  def process_output_commands(self,attrname='dataextraction'):
    """Process a list of data extraction commands
    Arguments:
      attrname = attribute name of self.modelparams containing the command list
    No return value."""
    for cmd in getattr(self.modelparams,attrname,[]):
      #Function name and arguments
      funcname, kwargs = cmd
      #Call it
      try:
        getattr(self,funcname)(**kwargs)
      except Exception as einst:
        print("Exception occured in %s for command: %s"%(attrname,str(cmd)), file=sys.stderr)
        raise einst
    return

  def create_output(self):
    """Process the data extraction commands for the completed solution
    Required attributes:
      diskwrite = boolean, True to write the output files.
        If false, the output will only be stored in the object itself
    No return value.
    Output files may be generated."""
    #Output location(s)
    if not osp.isdir(self.outdir):
      os.makedirs(self.outdir)

    #Process each command
    self.process_output_commands('dataextraction')

    #Put results into info
    self.info['results']=self.results
    #If present, add parametric locations
    if hasattr(self,'parametric_locations'):
      self.info['parametric_locations']=self.parametric_locations
    #Write output files if requested
    if self.diskwrite:
      common.writeyaml(self.info,osp.join(self.outdir,'info.yaml'))
      self.outdata.to_pickle(osp.join(self.outdir,'outdata.pkl'))

    #Done
    return

  ##TODO: support index argument for split fields
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

  def splitfield(self,namewhole,namesplit):
    """Call the split() method of a solution field.
    Arguments:
      namewhole = name of the field to split
      namesplit = attribute name to store the result in
    Required attributes:
      The attribute specified by namewhole.
    New attributes:
      The attribute specified by namesplit.
    No return value.
    No other side-effects."""
    setattr(self,namesplit,getattr(self,namewhole).split())

  def facet_area(self,pfacet,name,internal=False):
    """Compute the area of the specified facet.
    Arguments:
      pfacet = physical facet number of facet to calculate area of
      name = name for storage in the results dictionary
      internal = boolean, default False, True for internal boundary, False for external
    Required attributes:
      mesh = FEniCS Mesh object
      facet = FEniCS MeshFunction object for facet numbers
    No new attributes.
    New item(s) added to results dictionary.
    No return value."""
    if internal:
      integral_type='interior_facet'
    else:
      integral_type='exterior_facet'
    this_ds=fem.Measure(integral_type,domain=self.mesh,subdomain_data=self.facets)
    calcarea=fem.assemble(fem.Constant(1)*this_ds(pfacet))
    self.results[name]=calcarea
    return
