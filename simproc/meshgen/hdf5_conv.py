"""Convert XML mesh and mesh functions to HDF5 format"""

#This package
from ..requesthandler.filepath import Path
from ..requesthandler.request import Request
from ..requesthandler.yaml_manager import register_classes, read as readyaml
from ..requesthandler import locators

#Site packages
import fenics as fem

#Locators
locators.folder_structure.update(mesh_hdf5file=['mesh',0,'hdf5'])

#Constants
DCONV_XML_SUFFIX={'facet_xml':'_facet_region', 'cell_xml':'_physical_region'}

_HDF5Converter_props_schema_yaml="""#HDF5Converter
name: {type: string}
mesh_xml: {type: pathlike}
facet_xml: {type: pathlike}
cell_xml: {type: pathlike}
hdf5file: {type: pathlike}"""

class HDF5Converter(Request):
  """Convert FEniCS Mesh and MeshFunctions from XML to HDF5 format
  
  User-defined attributes:
  
    - mesh_xml = Path to input .xml file
        This is the path to output file containing the mesh itself.
    - facet_xml = optional Path to input .xml file containing facet meshfunction data
        Computed from the mesh_xml attribute if not provided, assuming the same directory
    - cell_xml = optional Path to input .xml file containing cell meshfunction data
        Computed from the mesh_xml attribute if not provided, assuming the same directory
    - hdf5file = Path to output .hdf5 file
  """
  _self_task=True
  _required_attrs=['mesh_xml','hdf5file']
  _outputfile_attrs=['hdf5file']
  _inputfile_attrs=['mesh_xml','facet_xml','cell_xml']
  _props_schema=readyaml(_HDF5Converter_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(HDF5Converter, self).__init__(**kwargs)
    #Compute paths for facet_xml and cell_xml if not provided
    mesh_xml=Path(self.mesh_xml)
    for attrname in DCONV_XML_SUFFIX.keys():
      if not hasattr(self,attrname):
        filename=mesh_xml.stemname+DCONV_XML_SUFFIX[attrname]+mesh_xml.ext
        setattr(self,attrname,Path(mesh_xml.folder,filename))
  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Read in the xml files (Mesh and MeshFunctions)
    mesh =   fem.Mesh(str(self.mesh_xml))
    facets = fem.MeshFunction("size_t", mesh, str(self.facet_xml))
    cells =  fem.MeshFunction("size_t", mesh, str(self.cell_xml))
    #Output to HDF5
    hdf5=fem.HDF5File(mesh.mpi_comm(),str(self.hdf5file),'w')
    hdf5.write(mesh,'mesh')
    hdf5.write(facets,'facets')
    hdf5.write(cells,'cells')
    hdf5.close()
    #Done
    return

#Register for loading from yaml
register_classes([HDF5Converter])
