"""Base functionality for simulation requests"""

#Standard library
import sys
from copy import deepcopy

#Site packages
import numpy as np
import fenics as fem

#This package
from ..requesthandler.customization import CustomizableRequest, make_schema
from ..requesthandler import yaml_manager
from ..requesthandler import locators
from .meshinfo import MeshInfo
from .plotseries import PlotSeries

#Locators
locators.folder_structure.update(SolutionFile=['solutions',0,1])

#To help with Conditions schemas

def update_schema_props(origschema,newprops,newreq=None):
  """Returns a new schema dictionary by adding to an existing one.
  
  Arguments:
  
    - origschema = the original schema (not just the properties), as a dictionary
    - newprops = the dictionary to be added to the properties key of the original schema
    - newreq = optional, sequence of properties to be added as required in the new schema
    
  Note that the original schema dictionary is not modified: a deep copy is created, modified, then returned."""
  newschema=deepcopy(origschema)
  newschema['properties'].update(newprops)
  if newreq is not None:
    newschema['required']+=newreq
  return newschema

def update_conditions(origschema,newconditions):
  """Update the conditions portion of a SimulationRequest schema"""
  newschema=deepcopy(origschema)
  newschema['conditions'].update(newconditions)
  return newschema

EmptyConditions_schema_yaml="""#EmptyConditions
type: object
properties: {}
required: []
additionalProperties: False
"""
EmptyConditions_schema=yaml_manager.readstring(EmptyConditions_schema_yaml)

GenericConditions_props_schema_yaml="""#GenericConditions
elementorder: {type: integer}
dirichlet: {type: object}
neumann: {type: object}
"""
GenericConditions_props_schema=yaml_manager.readstring(GenericConditions_props_schema_yaml)
GenericConditions_schema=update_schema_props(EmptyConditions_schema,GenericConditions_props_schema,['elementorder'])

_SimulationRequest_props_schema_yaml="""#SimulationRequest
mesh:
  anyOf:
    - {type: 'null'}
    - {type: pathlike}
meshmeta:
  anyOf:
    - {type: 'null'}
    - {type: pathlike}
hasmeshfuncs:
  type: boolean
conditions:
  type: object
dataextraction:
  type: array
loaddata:
  type: array
metadata:
  type: object
solver_parameters:
  type: object
meshinfo: {}
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
    - hasmeshfuncs: optional, False to avoid loading mesh functions
    - conditions: dictionary specifying model conditions such as element order, boundary conditions, etc.
    - dataextraction = a sequence of data extraction commands to execute after solving the model

      Each command is a pair (cmdname, arguments), where:

        - cmdname = name of the simulator object's method to call, as a string
        - arguments = dictionary of arguments to the method: {argname: value,...}

    - loaddata = a sequence of data loading commands to execute at simulator initialization
    
      Each command is a triple (attrpath, filepath, fieldtag), where:
      
        - attrpath = the attribute path of the simulator attribute, as a string, to store the loaded data
        - filepath = path to the data file to load, as string, relative to folderstructure.datafolder
        - fieldtag = string identifying the HDF5 field name to load
    
    - metadata = dictionary of metadata about the model, for use in post-processing"""
  _self_task=True
  _required_attrs=['name','mesh','conditions']
  _inputfile_attrs=['mesh','meshmeta']
  _config_attrs=['mesh','meshmeta','hasmeshfuncs','conditions','dataextraction','loaddata','metadata']
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
    self._more_outputfiles+=self.list_outputfiles(getattr(self,'dataextraction',[]))

  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Load the mesh, unless provided some other way
    if self.mesh is None:
      assert getattr(self,'meshinfo',None) is not None, "Must provide either MeshInfo or a mesh file to load."
      ##TODO: dependency checking on the mesh won't work for this, of course
    else:
      #Process mesh-related attributes
      meshmeta=getattr(self,'meshmeta',None)
      if meshmeta is not None:
        meshmeta=self.render(meshmeta)
      hasmeshfuncs=getattr(self,'hasmeshfuncs',True)
      #Load
      self.meshinfo=MeshInfo.load(self.render(self.mesh),meshmeta,hasmeshfuncs)
    #Do the simulation
    self.run_sim()
    #Generate output, if any
    self.process_command_sequence(attrpath='dataextraction',singlefunc=None,positional=False)
    return

  def run_sim(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'run_sim' method."%str(type(self)))

  def loadfield(self,attrpath,infpath,fieldtag,idx=None):
    """Load data into the simulator from an HDF5 input file
    
    Arguments:
    
      - attrpath = attribute path to load the data into, as string
      
        Note that this attribute must already exist, and be of the proper type to receive the requested data.
      
      - infpath = path to input file
      
      - fieldtag = string identifying the HDF5 field name to load
      
      - idx = index number (as integer) speciying location within the given attribute, None (default) if not the attribute itself is not a sequence
    
    No return value."""
    inputval = self.get_nested(attrpath)
    if idx is not None:
      inputval = inputval[idx]
    hdf5=fem.HDF5File(self.meshinfo.mesh.mpi_comm(),str(self.render(infpath)),'r')
    hdf5.read(inputval,fieldtag)
    hdf5.close()

  def process_load_commands(self):
    """Process a list of load data commands

    No return value."""
    self.process_command_sequence(attrpath='loaddata',singlefunc='loadfield',positional=True)
    return

  def list_outputfiles(self,cmdlist,filearg_list=None):
    """Get a list of all the files generated by the data extraction commands.

    Arguments:

      - cmdlist = list of data extraction commands,
        each command consists of pair (cmdname, arguments):

          - cmdname = name of data extraction method of the simulator class
          - arguments = dictionary of all arguments needed by the extraction method

    Return:

      - outfiles = list of generated output files (names only, not including their folder)"""
    #List of arguments that may contain output filenames
    if filearg_list is None:
      filearg_list=['filename','outfpath']
    #Initalize the output list
    outfiles=[]
    #Go through all commands in the list
    for cmdname, arguments in cmdlist:
      helperfunc=cmdname+'_outputfiles'
      if hasattr(self,helperfunc):
        outfiles += getattr(self,helperfunc)(**arguments)
      else:
        #Check all possible arguments that could contain the name of an output file
        present_args=[n for n in filearg_list if n in arguments.keys()]
        outfiles += [arguments[n] for n in present_args]
    return outfiles

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
    yaml_manager.writefile(outdict,self.render(outfpath))
    return

  def writefield_outputfiles(self,outfpath,attrpath='soln',idx=None,outname=None):
    """Compute a list of output files generated by writefield"""
    outlist=[outfpath]
    ofp=self.render(outfpath)
    if ofp.suffix.lower()=='.pvd':
      #There will also be .vtu or .pvtu files
      ##TODO: detect if MPI, and compute pvtu files
      vtufname=ofp.stemname + "%06d"%0 + ".vtu"
      vtupath=ofp.folder_path / vtufname
      outlist.append(vtupath)
    return outlist

  def writefield(self,outfpath,attrpath='soln',idx=None,outname=None):
    """Write field to VTK or HDF5 file

    Arguments:

      - outfpath = path to output file

        If the filename ends with ".hdf5" (case insensitive), an HDF5 file is created.
        Otherwise, the filetype is selected by FEniCS.

      - attrpath = attribute path to output, as string, defaults to 'soln'
      
      - idx = index number (as integer) of the solution field to write out, None (default) if not a sequence
      
      - outname = optional output field name within output file, as string, supported only for HDF5 format.
      
        If not provided, output field name defaults to value calculated from attrpath and idx.

    Required attributes:

      - the specified attribute containing the FEniCS field

    No new attributes.

    No return value.

    Output file is written."""
    output = self.get_nested(attrpath)
    if idx is not None:
      output = output[idx]
    if str(self.render(outfpath))[-5:].lower()=='.hdf5':
      #HDF5 format
      if outname is None:
        fieldtag=attrpath
        if idx is not None:
          fieldtag+='_%d'%idx
      else:
        fieldtag = outname
      hdf5=fem.HDF5File(self.meshinfo.mesh.mpi_comm(),self.renderstr(outfpath),'w')
      hdf5.write(output,fieldtag)
      hdf5.close()
    else:
      #Format controlled by FEniCS (including VTK files: .pvd, etc.)
      out_file=fem.File(str(self.render(outfpath)))
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

  def domain_volume(self,attrpath='volume',dxname='dx'):
    """Get the domain volume from integration
    
    Arguments:
    
      - attrpath = optional, attribute path for storing results, as string
      - dxname = optional, name of attribute with fenics domain volume measure, defaults to "dx"
    
    New attribute added/overwritten.
    No return value.
    No output files."""
    dx=getattr(self,dxname)
    volume=fem.assemble(fem.Constant(1)*dx)
    self.set_nested(attrpath,volume)

  def facet_area(self,pfacet,attrpath,internal=False):
    """Compute the area of the specified facet.

    Arguments:

      - pfacet = physical facet number of facet to calculate area of
      - attrpath = attribute path for storage of result
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

  def line_profile(self,startloc,endloc,num,plotpath,label,attrpath='soln',indep=None,idx=None):
    """Get data to plot a result along the specified line at a single point in time
  
    Arguments:
  
      - startloc = argument to get_pointcoords for start of line
      - endloc = argument to get_pointcoords for end of line
      - num = number of sampled points
      - indep = index of the coordinate parameter to use as the independent variable for the plot (zero-based) (omit to use distance from start point)
      - plotpath = attribute path for storing the generated PlotSeries instance
      - label = series label to assign, as string
      - attrpath = attribute path to data source, defaults to 'soln'
      - indep = identifier for independent variable:
          integer 0-d to use that coordinate of the point, or
          None (default) to use distance from the start point
      - idx = index of the solution field to write out, None (default) if not a sequence
  
    No return value."""
    #Get the object with the data
    vals=self.get_nested(attrpath)
    if idx is not None:
      vals = vals[idx]
    #Get the points for data extraction
    assert len(startloc)==len(endloc), "Start and end locations have different dimensionality"
    startcoords=self.get_pointcoords(startloc)
    endcoords=self.get_pointcoords(endloc)
    start_ends=[itm for itm in zip(startcoords,endcoords)]
    ranges=[np.linspace(start,end,num) for start,end in start_ends]
    points=[t for t in zip(*ranges)]
    #Function to calculate independent variable for a given point
    if indep is None:
      indep_f = lambda pt: np.sqrt(sum([(startcoords[i]-c)**2 for i,c in enumerate(pt)]))
    else:
      indep_f = lambda pt: pt[indep]
    #Extract data points
    llist=[]
    vlist=[]
    for pt in points:
      try:
        vlist.append(vals(*pt))
        llist.append(indep_f(pt))
      except RuntimeError:
        pass #point is not inside mesh; skip
    #Create PlotSeries
    larr=np.array(llist)
    varr=np.array(vlist)
    series=PlotSeries(xvals=larr,yvals=varr,label=label)
    #Store data
    self.set_nested(plotpath,series)
    return

#Register for loading from yaml
yaml_manager.register_classes([SimulationRequest])

