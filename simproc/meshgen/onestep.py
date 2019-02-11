"""Run gmsh, dolfin-convert, and the HDF5 conversion with a single request"""

#This package
from . import gmsh_runner
from . import dconv_runner
from . import hdf5_conv


_GeoToHDF5Request_props_schema_yaml="""#GeoToHDF5Request
name: {type: string}
geofile: {type: pathlike}
mshfile: {type: pathlike}
gmsh_txtfile: {type: pathlike}
meshmetafile: {type: pathlike}
mshfile: {type: pathlike}
dconv_txtfile: {type: pathlike}
mesh_xml: {type: pathlike}
facet_xml: {type: pathlike}
cell_xml: {type: pathlike}
hdf5file: {type: pathlike}"""

class GeoToHDF5Request(Request):
  """Run gmsh, dolfin-convert, and the HDF5 conversion with a single request
  
  User-defined attributes:

    - geofile: Path to input .geo file
    - mshfile: Path to output .msh file
    - gmsh_txtfile: Path to text file to store gmsh message output
    - meshmetafile: optional, Path to yaml file to store mesh metadata
    - integer_arg: optional, integer to pass to gmsh on the command line, to specify meshing dimension
        defaults to 0, which indicates that the .geo file contains the appropriate `Mesh` command.
    - dconv_txtfile: Path to text file to store dolfin-convert message output
    - mesh_xml = Path to .xml file for the mesh itself
    - facet_xml = Path to .xml file containing facet meshfunction data
    - cell_xml = Path to .xml file containing cell meshfunction data
    - hdf5file = Path to output .hdf5 file

  Calculated Attributes:
  
    - gmsh_request: Request that will run gmsh
    - dconv_request: Request that will run dolfin-convert
    - hdf5_request: Request that will convert to hdf5"""
  _required_attrs=['geofile']
  _outputfile_attrs=['mshfile','gmsh_txtfile','meshmetafile','dconv_txtfile','mesh_xml','facet_xml','cell_xml','hdf5file']
  _inputfile_attrs=['geofile']
  _child_attrs=['gmsh_request','dconv_request','hdf5_request']
  _props_schema=readyaml(_GeoToHDF5Request_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(GeoToHDF5Request, self).__init__(**kwargs)
    #Set up child requests
    self.gmsh_request=gmsh_runner.GmshRequest(name=self.name+".gmsh",
                                              geofile=self.geofile,
                                              mshfile=self.mshfile,
                                              txtfile=self.gmsh_txtfile,
                                              meshmetafile=self.meshmetafile,
                                              integer_arg=self.integer_arg)
    self.dconv_request=dconv_runner.DolfinConvertRequest(name=self.name+'.dconv',
                                                         mshfile=self.mshfile,xmlfile=self.mesh_xml,txtfile=self.dconv_txtfile)
    self.hdf5_request=hdf5_conv.HDF5ConvertRequest(name=self.name+'.hdf5',
                                                   mesh_xml=self.mesh_xml,facet_xml=self.facet_xml,cell_xml=self.cell_xml,
                                                   hdf5file=self.hdf5file)

#Register for loading from yaml
register_classes([GeoToHDF5Request])

