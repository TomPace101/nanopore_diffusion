"""Base functionality for simulation requests"""

#Site packages
import fenics as fem

#This package
from ..requesthandler.customization import CustomizableRequest
from .meshinfo import MeshInfo

_SimulationRequest_props_schema_yaml="""#SimulationRequest
mesh:
  anyOf:
    - type: pathlike
    - type: object
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
  
  User-defined attributes:
  
    - mesh: path to mesh hdf5 file, or dictionary specifying mesh, depending on the simulation
    - conditions: dictionary specifying model conditions such as element order, boundary conditions, etc.
    - dataextraction: TODO
    - loaddata: TODO
    - metadata: TODO"""##TODO
  _self_task=True
  _required_attrs=[] ##TODO
  _outputfile_attrs=[] ##TODO
  _inputfile_attrs=[] ##TODO
  _props_schema=readyaml(_SimulationRequest_props_schema_yaml)

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

  def facet_area(self,pfacet,name,internal=False):
    """Compute the area of the specified facet.

    Arguments:

      - pfacet = physical facet number of facet to calculate area of
      - name = name for storage in the results dictionary
      - internal = boolean, default False, True for internal boundary, False for external

    Required attributes:

      - mesh = FEniCS Mesh object
      - facet = FEniCS MeshFunction object for facet numbers

    No new attributes.

    New item(s) added to results dictionary.

    No return value."""
    if internal:
      integral_type='interior_facet'
    else:
      integral_type='exterior_facet'
    this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
    calcarea=fem.assemble(fem.Constant(1)*this_ds(pfacet))
    self.results[name]=calcarea
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

  def line_profile(self,startloc,endloc,num,plotname,label,attrname='soln',indep=None,idx=None):
    """Get data to plot a result along the specified line at a single point in time

    Arguments:

      - startloc = argument to get_pointcoords for start of line
      - endloc = argument to get_pointcoords for end of line
      - num = number of sampled points
      - indep = index of the coordinate parameter to use as the independent variable for the plot (zero-based) (omit to use distance from start point)
      - plotname = name of plot in outdata.plots, as string
      - label = series label to assign, as string
      - attrname = name of attribute to output, as string, defaults to 'soln'
      - indep = identifier for independent variable:
          integer 0-d to use that coordinate of the point, or
          None (default) to use distance from the start point
      - idx = index of the solution field to write out, None (default) if not a sequence

    Required attributes:

      - outdata = instance of OutData
      - mesh_metadata = only required if needed by location specifiers, dictionary of mesh metadata

    No new attributes.

    Nothing added to results dictionary.

    No return value.

    Series is added to ``outdata.plots``."""
    #Get the object with the data
    vals=getattr(self,attrname)
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
    series=plotdata.PlotSeries(xvals=larr,yvals=varr,label=label)
    #Store data
    if not plotname in self.outdata.plots.keys():
      self.outdata.plots[plotname]=[]
    self.outdata.plots[plotname].append(series)
    return

#Register for loading from yaml
register_classes([SimulationRequest])

