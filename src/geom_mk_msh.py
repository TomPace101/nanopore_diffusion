#Call gmsh to convert .geo file(s) into .msh, based on derivative of MeshParameters

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp
from subprocess import call
import sys

#Site packages

#Local
import folderstructure as FS
import common

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

#Template for gmsh command
#arguments: mshfile, geofile, txtfile
cmd_tmpl="gmsh -0 -setstring paramlocsfile %s -o %s %s >%s"

class GmshRunner(common.ParameterSet):
  __slots__=('meshname','geomdefname','tmplvalues','geofile','mshfile','txtfile','paramlocsfile','_folders')
  _required_attrs=['meshname','geomdefname','tmplvalues']
  _config_attrs=_required_attrs
  #don't need sourcefile as input file due to config
  _inputfile_attrs=['geofile']
  _more_inputfiles=[thisfile,common.__file__]
  _outputfile_attrs=['mshfile','txtfile','paramlocsfile']
  _taskname_src_attr='meshname'
  #Remember loaded geometry definitions so we aren't reading the same files over and over when generating multiple meshes
  loaded_geomdefs={}

  def __init__(self,**kwd):
    #Initialization from base class
    super(GmshRunner, self).__init__(**kwd)
    #Get folders
    self._folders={'geofile':osp.join(FS.geofolder,self.basename),
                   'mshfile':osp.join(FS.mshfolder,self.basename),
                   'txtfile':osp.join(FS.gmsh_outfolder,self.basename),
                   'paramlocsfile':osp.join(FS.paramlocs_outfolder,self.basename)} 
    #Get name of input and output files
    self.geofile=self.meshname+'.geo'
    self.mshfile=self.meshname+'.msh'
    self.txtfile=self.meshname+'.txt'
    self.paramlocsfile=self.meshname+'.yaml'

  def run(self):
    print(self.mshfile)
    #Create directories if necessary
    for oattr in self._outputfile_attrs:
      if not osp.isdir(self._folders[oattr]):
        os.makedirs(self._folders[oattr])
    #Run the shell command
    cmd_str=cmd_tmpl%(self.full_path('paramlocsfile'),self.full_path('mshfile'),self.full_path('geofile'),self.full_path('txtfile'))
    call(cmd_str,shell=True)

#Support command-line arguments
if __name__ == '__main__':
  program_description='Create gmsh .msh file(s) from .geo file(s) by running gmsh from a yaml input file'
  input_file_description="""Path to parameter definition file for the mesh
    This is a potentially multi-doc yaml file, where each document specifies one mesh to generate."""
  
  common.run_cmd_line(program_description,input_file_description,GmshRunner)
