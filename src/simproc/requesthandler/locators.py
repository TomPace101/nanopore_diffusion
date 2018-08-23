"""Define data locators

This module has a variable ``TOPFOLDER`` which specifies the path of the folder to work in.
All locators return paths that start at this folder.

From within python, you can simply set ``TOPFOLDER`` as a variable,
and all locators will return paths based on the new location.

For scripts, use the environment variable ``TOPFOLDER`` to tell the program where this folder is."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
import os

#Site packages

#This package
from . import filepath
from . import yaml_manager

#Constants
topfolder_environ='TOPFOLDER'

#Get path to the top folder
if topfolder_environ in os.environ.keys():
  # TOPFOLDER=Path(osp.normpath(osp.abspath(os.environ[topfolder_environ])))
  TOPFOLDER=filepath.Path(os.environ[topfolder_environ]).expanduser().resolve()
else:
  modpath=filepath.Path(__file__)
  srcfolder=modpath.parent.parent.parent
  TOPFOLDER=srcfolder.parent / 'data'

class DataFile(object):
  """File location relative to ``TOPFOLDER``
  
  The location of ``TOPFOLDER`` is specified by the environment variable ``TOPFOLDER``.
  If this environment variable does not exist, a default value is provided.
  See ``folderstructure.py``."""
  def __init__(self,*args,**kwargs):
    self.subpath=filepath.Path(*args,**kwargs)
  def path(self,req):
    return TOPFOLDER / self.subpath
  @classmethod
  def from_yaml(cls, constructor, node):
    return cls(node.value)

def locator_factory(ltype):
  """Factory function to return a locator class for a given name
  
  The locator class returned will construct file paths as follows:
  
    - all paths begin with ``locators.TOPFOLDER``
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
    def path(self,reqname):
      namelist=reqname.split('.')
      specifier=folder_structure[ltype]
      out=TOPFOLDER
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
  Hence, this is an object that simply calls another function when initialized,
  and then does nothing else ever."""
  def __init__(self,**kwargs):
    folder_structure.update(kwargs)

#Register for loading from yaml
yaml_manager.register_classes([DataFile,UpdateFolderStructure])
