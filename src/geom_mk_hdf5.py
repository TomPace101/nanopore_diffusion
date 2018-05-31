"""Convert XML mesh and mesh functions to HDF5 format"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp
import sys

#Site packages
import fenics as fem

#Local
import folderstructure as FS
import common

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

class HDF5Converter(common.ParameterSet):
  """Subclass of common.ParameterSet for converting FEniCS Mesh and MeshFunctions from XML to HDF5 format
  
  Attributes:
  
    To be read in:
    
      - meshname, geomdefname, tmplvalues as for buildgeom.MeshParameters
    
    To be generated by methods:
    
      - xmlfolder = path to xml input files for mesh
      - mesh_xml, facet_xml, cell_xml = various mesh input xml files, full paths, as strings
      - hdf5file = file name for hdf5 file (not full path)
  """
  __slots__=('meshname','geomdefname','tmplvalues','_folders','xmlfolder','mesh_xml','facet_xml','cell_xml','hdf5file')
  _required_attrs=['meshname','geomdefname','tmplvalues']
  _config_attrs=_required_attrs
  #don't need sourcefile as input file due to config
  _inputfile_attrs=['mesh_xml', 'facet_xml', 'cell_xml']
  _more_inputfiles=[thisfile,common.__file__]
  _outputfile_attrs=['hdf5file']
  _taskname_src_attr='meshname'

  def __init__(self,**kwd):
    #Initialization from base class
    super(HDF5Converter, self).__init__(**kwd)
    #Get folders
    self._folders={'hdf5file':osp.join(FS.mesh_hdf5_folder,self.basename)}
    #Input files
    self.xmlfolder=osp.join(FS.xmlfolder,self.basename)
    self.mesh_xml=osp.join(self.xmlfolder,self.meshname+'.xml')
    self.facet_xml=osp.join(self.xmlfolder,self.meshname+'_facet_region.xml')
    self.cell_xml=osp.join(self.xmlfolder,self.meshname+'_physical_region.xml')
    #Output files
    self.hdf5file=self.meshname+'.hdf5'
  
  def run(self):
    print(self.hdf5file)
    #Create directories if necessary
    for oattr in self._outputfile_attrs:
      if not osp.isdir(self._folders[oattr]):
        os.makedirs(self._folders[oattr])
    #Read in the xml files (Mesh and MeshFunctions)
    mesh =  fem.Mesh(self.mesh_xml)
    facets =fem.MeshFunction("size_t", mesh, self.facet_xml)
    cells = fem.MeshFunction("size_t", mesh, self.cell_xml)
    #Output to HDF5
    hdf5=fem.HDF5File(mesh.mpi_comm(),self.full_path('hdf5file'),'w')
    hdf5.write(mesh,'mesh')
    hdf5.write(facets,'facets')
    hdf5.write(cells,'cells')

#Support command-line arguments
if __name__ == '__main__':
  program_description='Create dolfin mesh .hdf5 file(s) from .xml file(s)'
  input_file_description="""Path to parameter definition file for the mesh
    This is a potentially multi-doc yaml file, where each document specifies one mesh to generate."""
  
  common.run_cmd_line(program_description,input_file_description,HDF5Converter)