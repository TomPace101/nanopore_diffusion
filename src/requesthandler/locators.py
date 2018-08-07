"""Define data locators"""
##TODO: a way to redefine datafolder from a yaml file

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
import os

#Site packages

#This package
from . import filepath
from . import yaml_manager

#Locate data folder
if 'DATALOC' in os.environ.keys():
  # datafolder=Path(osp.normpath(osp.abspath(os.environ['DATALOC'])))
  datafolder=filepath.Path(os.environ['DATALOC']).expanduser().reslolve()
else:
  srcfolder=filepath.Path(__file__).parent.parent
  datafolder=srcfolder.parent / 'data'


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
  """Factory function to return a locator class for a given name
  
  The locator class returned will construct file paths as follows:
  
    - all paths begin with ``locators.datafolder``
    - then, one subdirectory is added for each element of the folder structure definition sequence for the locator (see below)
    - finally, the subpath used to initialize the locator is appended
  
  The folder structure definition sequence defines subdirectories as follows:
  
    - any element which is a string adds the string as a folder name
    - any element which is an integer adds the matching element in the dotted name of the request
    - integers beyond the last element of the dotted name are simply ignored.
  
  Arguments:
  
    - ltype = name of locator class, as string
  
  Returns:
  
    - lclass = locator class"""
  class lclass(DataFile):
    """File location as specified by the name of its parent request and the current folder structure settings"""
    def __init__(self,*args,**kwargs):
      self.subpath=filepath.Path(*args,**kwargs)
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
      out /= self.subpath
      return out
  lclass.__name__=ltype
  return lclass

class FolderStructure(odict):
  """For storing an expected folder structure.
  
  This dictionary is intended to have the following organization: {locator_class_name: [structure_defintion_sequence]}  
  
    - each key is a string, matching the name of a locator class
    - the structure definition sequence indicates how to create the path for files using the locator class
    
  See the locator class factory for explanation of how the structure definition sequence determines the file path.""" 
  def update(self,**kwargs):
    #Before we actually change the folder structure, we need to check for new locator classes
    #All existing locator types
    existing=self.keys()
    #Newly defined types
    newtypenames=[lt for lt in kwargs.keys() if lt not in existing]
    newtypes=[]
    for ltype in newtypenames:
      #Create the locator
      lclass=locator_factory(ltype)
      globals()[ltype]=lclass
      newtypes.append(lclass)
    #Register the new locators for reading from yaml
    yaml_manager.register_classes(newtypes)
    #Update the expected folder structure
    super(FolderStructure, self).update(**kwargs)

#Folder structure singleton
folder_structure=FolderStructure()
# folder_structure={
#   'CustomizationFile':['customizations'],
#   'RequestFile':['requests'],
#   'MeshGeomdefFile':['mesh','geomdef'],
#   'MeshTemplateFile':['mesh','templates'],
#   'MeshGeoFile':[0,'mesh',1],
#   'MeshMshFile':[0,'mesh',1],
#   'MeshXmlFile':[0,'mesh',1],
#   'MeshHdf5File':[0,'mesh',1],
#   'MeshGmshOutFile':[0,'mesh',1],
#   'MeshMetadataFile':[0,'mesh',1],
#   'SolutionFile':[0,1],
#   'PostprocFile':[0,1]
# }

class UpdateFolderStructure(object):
  """Apply changes to the folder structure used to locate files
  
  This isn't really a class. It's a function.
  But when loading from yaml, there isn't a way to call a funtion directly.
  You can only load classes.
  Hence, this hack of an object that simply calls another function when initialized,
  and then does nothing else ever."""
  def __init__(self,**kwargs):
    folder_structure.update(kwargs)

#Register for loading from yaml
yaml_manager.register_classes([DataFile,UpdateFolderStructure])
