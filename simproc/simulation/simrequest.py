"""Base functionality for simulation requests"""

#Site packages
import fenics as fem

#This package
from ..requesthandler.customization import CustomizableRequest, make_schema
from ..requesthandler import yaml_manager
from .meshinfo import MeshInfo

def list_outputfiles(cmdlist):
  """Get a list of all the files generated by the data extraction commands.

  Arguments:

    - cmdlist = list of data extraction commands,
      each command consists of pair (cmdname, arguments):

        - cmdname = name of data extraction method of the simulator class
        - arguments = dictionary of all arguments needed by the extraction method

  Return:

    - outfiles = list of generated output files (names only, not including their folder)"""
  #This is a bit of a guess: we assume we know the argument names that can hold output file paths
  filearg_list=['filename','outfpath']
  outfiles=[]
  for cmdname, arguments in cmdlist:
    #Check all possible arguments that could contain the name of an output file
    present_args=[n for n in filearg_list if n in arguments.keys()]
    outfiles.extend([arguments[n] for n in present_args])
  return outfiles

_SimulationRequest_props_schema_yaml="""#SimulationRequest
mesh:
  type: pathlike
meshmeta:
  type: pathlike
conditions:
  type: object
dataextraction:
  type: array
loaddata:
  type: array
metadata:
  type: object
"""

class SimulationRequest(CustomizableRequest):
  """Base class for FEniCS simulations
  
  The `run` method of this class includes only basic functionality:

    - Do the standard pre-run checks
    - Load the mesh
    - Call the method `run_sim`, which is to be defined by derived classes or customization
    - Process the output commands.
  
  Of course, even this simple behavior can be overidden through customization if desired.
  
  User-defined attributes:
  
    - mesh: path to mesh hdf5 file
    - meshmeta: optional, path to mesh metadata yaml file
    - conditions: dictionary specifying model conditions such as element order, boundary conditions, etc.
    - dataextraction = a sequence of data extraction commands to execute after solving the model

      Each command is a pair (cmdname, arguments), where:

        - cmdname = name of the simulator object's method to call, as a string
        - arguments = dictionary of arguments to the method: {argname: value,...}

    - loaddata = a sequence of data loading commands to execute at simulator initialization
    
      Each command is a triple (attrname, filepath, fieldtag), where:
      
        - attrname = the name of the simulator attribute, as a string, to store the loaded data
        - filepath = path to the data file to load, as string, relative to folderstructure.datafolder
        - fieldtag = string identifying the HDF5 field name to load
    
    - metadata = dictionary of metadata about the model, for use in post-processing"""
  _self_task=True
  _required_attrs=['mesh','conditions']
  _inputfile_attrs=['mesh','meshmeta']
  _props_schema=make_schema(_SimulationRequest_props_schema_yaml)

  def __init__(self,**kwargs):
    #Initialization from base class
    super(SimulationRequest, self).__init__(**kwargs)
    #Get input files from load data commands
    self._more_inputfiles=getattr(self,'_more_inputfiles',[]) #Initialize attribute if it doesn't already exist
    for aname,fpath,ftag in getattr(self,'loaddata',[]):
      self._more_inputfiles.append(fpath)
    #Get output files from data extraction commands
    self._more_outputfiles=getattr(self,'_more_outputfiles',[]) #Initialize attribute if it doesn't already exist
    self._more_outputfiles+=list_outputfiles(getattr(self,'dataextraction',[]))

  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Load the mesh
    meshmeta=getattr(self,'meshmeta',None)
    self.meshinfo=Meshinfo(self.mesh,meshmeta)
    #Do the simulation
    self.run_sim()
    #Generate output, if any
    self.process_output_command('dataextraction')
    return

  def run_sim(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'run_sim' method."%str(type(self)))

  def get_nested(self,dpath):
    """Return the value from the specified attribute/key/index path
    
    Arguments:
    
      - dpath = string describing path to the data, using dot separators, or a sequence
          The path may contain attributes and dictionary keys, with no need to distinguish between them.
          List indices are also allowed.
    
    Returns the requested data."""
    nxt=self
    if isinstance(dpath,str):
      seq=dpath.split('.')
    else:
      seq=dpath
    for name in seq:
      if hasattr(nxt,name):
        nxt = getattr(nxt,name)
      else:
        try:
          nxt.__getitem__(name)
        except:
          raise KeyError('Invalid path %s: No attribute, key, or index %s'%(dpath,name))
    return nxt

  def set_nested(self,dpath,val):
    """Set the value at the specified attribute/key/index path
    
    Arguments:
    
      - dpath = string describing path to the data, using dot separators, or a sequence
      - val = value to assign
    
    No return value."""
    if isinstance(dpath,str):
      seq=dpath.split('.')
    else:
      seq=dpath
    head=seq[:-1]
    tail=seq[-1]
    parent=self.get_nested(head)
    if hasattr(parent,'__setitem__'):
      parent.__setitem__(tail,val)
    else:
      setattr(parent,tail,val)
    return

  def loadfield(self,infpath,fieldtag,attrname,idx=None):
    """Load data into the simulator from an HDF5 input file
    
    Arguments:
    
      - infpath = path to input file, as string, relative to folderstructure.datafolder
      
      - fieldtag = string identifying the HDF5 field name to load
      
      - attrname = attribute name to load the data into, as string
      
        Note that this attribute must already exist, and be of the proper type to receive the requested data.
      
      - idx = index number (as integer) speciying location within the given attribute, None (default) if not the attribute itself is not a sequence
    
    No return value."""
    fullpath=osp.join(FS.datafolder,infpath)
    inputval = getattr(self,attrname)
    if idx is not None:
      inputval = inputval[idx]
    hdf5=fem.HDF5File(self.meshinfo.mesh.mpi_comm(),fullpath,'r')
    hdf5.read(inputval,fieldtag)
    hdf5.close()

  def process_load_commands(self,attrname='loaddata'):
    """Process a list of load data commands

    Arguments:

      - attrname = attribute name of self.modelparams containing the command list

    No return value."""
    for cmd in getattr(self.modelparams,attrname,[]):
      #Function arguments
      aname, fpath, ftag = cmd
      #Call it
      try:
        self.loadfield(fpath,ftag,aname)
      except Exception as einst:
        print("Exception occured in %s for command: %s"%(attrname,str(cmd)), file=sys.stderr)
        raise einst
    return

  def process_output_commands(self,attrname='dataextraction'):
    """Process a list of data extraction commands

    Arguments:

      - attrname = attribute name of self.modelparams containing the command list

    No return value."""
    for cmd in getattr(self,attrname,[]):
      #Function name and arguments
      funcname, kwargs = cmd
      #Call it
      try:
        getattr(self,funcname)(**kwargs)
      except Exception as einst:
        print("Exception occured in %s for command: %s"%(attrname,str(cmd)), file=sys.stderr)
        raise einst
    return

  def reportvalues(self,outfpath,mapping):
    """Write the selected output results to a yaml file
    
    Arguments:
    
      - outfpath = path to the output yaml file
      - mapping = mapping of output field names to object paths suitable for get_nested
    
    No attributes modified.
    Output file is created/overwritten.
    No return value."""
    outdict={}
    for key,dpath in mapping.items():
      outdict[key]=self.get_nested(dpath)
    yaml_manager.writefile(outdict,outfpath)
    return

  def writefield(self,outfpath,attrname='soln',idx=None,outname=None):
    """Write field to VTK or HDF5 file

    Arguments:

      - outfpath = path to output file

        If the filename ends with ".hdf5" (case insensitive), an HDF5 file is created.
        Otherwise, the filetype is selected by FEniCS.

      - attrname = name of attribute to output, as string, defaults to 'soln'
      
      - idx = index number (as integer) of the solution field to write out, None (default) if not a sequence
      
      - outname = optional output field name within output file, as string, supported only for HDF5 format.
      
        If not provided, output field name defaults to value calculated from attrname and idx.

    Required attributes:

      - the specified attribute containing the FEniCS field

    No new attributes.

    No return value.

    Output file is written."""
    output = getattr(self,attrname)
    if idx is not None:
      output = output[idx]
    if outfpath.suffix.lower()=='.hdf5':
      #HDF5 format
      if outname is None:
        fieldtag=attrname
        if idx is not None:
          fieldtag+='_%d'%idx
      else:
        fieldtag = outname
      hdf5=fem.HDF5File(self.meshinfo.mesh.mpi_comm(),outfpath,'w')
      hdf5.write(output,fieldtag)
      hdf5.close()
    else:
      #Format controlled by FEniCS (including VTK files: .pvd, etc.)
      out_file=fem.File(outfpath)
      out_file << output
    return

  def splitfield(self,namewhole,namesplit):
    """Call the split() method of a solution field.

    Arguments:

      - namewhole = name of the field to split
      - namesplit = attribute name to store the result in

    Required attributes:

      - The attribute specified by namewhole.

    New attributes:

      - The attribute specified by namesplit.

    No return value.

    No other side-effects."""
    setattr(self,namesplit,getattr(self,namewhole).split())

  def domain_volume(self,attrname='volume',dxname='dx'):
    """Get the domain volume from integration
    
    Arguments:
    
      - attrname = optional, name of attribute for storing results, as string
      - dxname = optional, name of attribute with fenics domain volume measure, defaults to "dx"
    
    New attribute added/overwritten.
    No return value.
    No output files."""
    dx=getattr(self,dxname)
    volume=fem.assemble(fem.Constant(1)*dx)
    setattr(self,attrname,volume)

  def facet_area(self,pfacet,attrname,internal=False):
    """Compute the area of the specified facet.

    Arguments:

      - pfacet = physical facet number of facet to calculate area of
      - attrname = attribute name for storage of result
      - internal = boolean, default False, True for internal boundary, False for external

    Required attributes:

      - mesh = FEniCS Mesh object
      - facet = FEniCS MeshFunction object for facet numbers

    New attribute added/overwritten.

    No return value."""
    if internal:
      integral_type='interior_facet'
    else:
      integral_type='exterior_facet'
    this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    calcarea=fem.assemble(fem.Constant(1)*this_ds(pfacet))
    setattr(self,name,calcarea)
    return

  def get_pointcoords(self,location):
    """Process a location specifier.

    Arguments:

      - location = specifier of location within the mesh

        This should be a tuple, with length matching the problem dimensions.
        Each entry is either a number or a string.
        Numbers represent physical coordinates within the mesh.
        Strings are replaced with the corresponding entry from the mesh metadata dictionary.
          (which is a required attributed in that case)

    Required attributes:

      - mesh_metadata = only required if needed by location specifiers, dictionary of mesh metadata

    Returns:

      - coords = the converted tuple"""
    coords=tuple()
    for v in location:
      if type(v)==str:
        v=self.meshinfo.metadata[v]
      coords+=(v,)
    return coords

  # def line_profile(self,startloc,endloc,num,plotname,label,attrname='soln',indep=None,idx=None):
  #   """Get data to plot a result along the specified line at a single point in time
  # 
  #   Arguments:
  # 
  #     - startloc = argument to get_pointcoords for start of line
  #     - endloc = argument to get_pointcoords for end of line
  #     - num = number of sampled points
  #     - indep = index of the coordinate parameter to use as the independent variable for the plot (zero-based) (omit to use distance from start point)
  #     - plotname = name of plot in outdata.plots, as string
  #     - label = series label to assign, as string
  #     - attrname = name of attribute to output, as string, defaults to 'soln'
  #     - indep = identifier for independent variable:
  #         integer 0-d to use that coordinate of the point, or
  #         None (default) to use distance from the start point
  #     - idx = index of the solution field to write out, None (default) if not a sequence
  # 
  #   Required attributes:
  # 
  #     - outdata = instance of OutData
  #     - mesh_metadata = only required if needed by location specifiers, dictionary of mesh metadata
  # 
  #   No new attributes.
  # 
  #   No return value.
  # 
  #   Series is added to ``outdata.plots``.""" ##TODO: we don't have outdata now
  #   #Get the object with the data
  #   vals=getattr(self,attrname)
  #   if idx is not None:
  #     vals = vals[idx]
  #   #Get the points for data extraction
  #   assert len(startloc)==len(endloc), "Start and end locations have different dimensionality"
  #   startcoords=self.get_pointcoords(startloc)
  #   endcoords=self.get_pointcoords(endloc)
  #   start_ends=[itm for itm in zip(startcoords,endcoords)]
  #   ranges=[np.linspace(start,end,num) for start,end in start_ends]
  #   points=[t for t in zip(*ranges)]
  #   #Function to calculate independent variable for a given point
  #   if indep is None:
  #     indep_f = lambda pt: np.sqrt(sum([(startcoords[i]-c)**2 for i,c in enumerate(pt)]))
  #   else:
  #     indep_f = lambda pt: pt[indep]
  #   #Extract data points
  #   llist=[]
  #   vlist=[]
  #   for pt in points:
  #     try:
  #       vlist.append(vals(*pt))
  #       llist.append(indep_f(pt))
  #     except RuntimeError:
  #       pass #point is not inside mesh; skip
  #   #Create PlotSeries
  #   larr=np.array(llist)
  #   varr=np.array(vlist)
  #   series=plotdata.PlotSeries(xvals=larr,yvals=varr,label=label) ##TODO
  #   #Store data
  #   if not plotname in self.outdata.plots.keys(): ##TODO we don't have outdata now
  #     self.outdata.plots[plotname]=[]
  #   self.outdata.plots[plotname].append(series)
  #   return

#Register for loading from yaml
yaml_manager.register_classes([SimulationRequest])

