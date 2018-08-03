"""Define data locators"""
##TODO: explain how locators work with folderstructure
##TODO: a way to redefine datafolder from a yaml file

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os

#Site packages

#This package
from . import filepath
from . import requestfile

#Locate data folder
if 'DATALOC' in os.environ.keys():
  # datafolder=Path(osp.normpath(osp.abspath(os.environ['DATALOC'])))
  datafolder=filepath.Path(os.environ['DATALOC']).expanduser().reslolve()
else:
  srcfolder=filepath.Path(__file__).parent.parent
  datafolder=srcfolder.parent / 'data'

#Default folder structure
folder_structure={
  'CustomizationFile':['customizations'],
  'RequestFile':['requests'],
  'MeshGeomdefFile':['mesh','geomdef'],
  'MeshTemplateFile':['mesh','templates'],
  'MeshGeoFile':[0,'mesh',1],
  'MeshMshFile':[0,'mesh',1],
  'MeshXmlFile':[0,'mesh',1],
  'MeshHdf5File':[0,'mesh',1],
  'MeshGmshOutFile':[0,'mesh',1],
  'MeshMetadataFile':[0,'mesh',1],
  'SolutionFile':[0,1],
  'PostprocFile':[0,1]
}

class DataFile(object):
  """File location relative to ``datafolder``
  
  The location of ``datafolder`` is specified by the environment variable ``DATALOC``.
  If this environment variable does not exist, a default value is provided.
  See ``folderstructure.py``."""
  parentpath=datafolder
  def __init__(self,*args,**kwargs):
    self.subpath=filepath.Path(*args,**kwargs)
  def path(self,req):
    return self.parentpath / self.subpath
  @classmethod
  def from_yaml(cls, constructor, node):
    return cls(node.value)

def locator_factory(ltype):
  class lclass(DataFile):
    """File location as specified by the name of its parent request and the current folder structure settings"""
    def __init__(self,*args,**kwargs):
      self.subpath=filepath(*args,**kwargs)
    def path(self,req):
      namelist=req.name.split('.')
      specifier=folder_structure[ltype]
      out=self.parentpath
      for itm in specifier:
        if type(itm) is int:
          if itm < len(namelist):
            out /= namelist[itm]
        else:
          out /= itm
      return out
  lclass.__name__=ltype
  return lclass

fs_locators=[]
for ltype in folder_structure.keys():
  lclass=locator_factory(ltype)
  globals()[ltype]=lclass
  fs_locators.append(lclass)

class UpdateFolderStructure(object):
  """Apply changes to the folder structure used to locate files
  
  The arguments to the class constructor are used to update the folder structure dictionary.
  This is a dictionary with keys matching the names of locator classes,
  and values specifying how to construct the expected file locations.
  
  For each key, provide a sequence of path components.
  folderstructure.datafolder is always prepended to the path.
  Each element of the sequence adds a subdirectory.
  Strings are added directly as subdirectories.
  Integers add subdirectories based on the request name.
  The request name is split by any dots,
  and the integers represent items in that sequence.
  Integers beyond the last sequence entry are simply ignored."""
  def __init__(self,**kwargs):
    folder_structure.update(kwargs)
    ##TODO: what if a new filetype was just define? We need to create a locator for it.
    ##And then that locator class needs to be registered for yaml loading too

#Register for loading from yaml
requestfile.register_classes([DataFile,UpdateFolderStructure]+fs_locators)
