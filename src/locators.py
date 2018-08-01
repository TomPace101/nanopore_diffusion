"""Define data locators"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#Local
import filepath
import folderstructure as FS

class DataFile(object):
  """File location relative to ``datafolder``
  
  The location of ``datafolder`` is specified by the environment variable ``DATALOC``.
  If this environment variable does not exist, a default value is provided.
  See ``folderstructure.py``."""
  parentpath=FS.datafolder
  def __init__(self,*args,**kwargs):
    self.subpath=filepath.Path(*args,**kwargs)
  def path(self,req):
    return self.parentpath / self.subpath

def locator_factory(ltype):
  class lclass(DataFile):
    """File location as specified by the name of its parent request and the current folder structure settings"""
    def __init__(self,*args,**kwargs):
      self.subpath=filepath(*args,**kwargs)
    def path(self,req):
      namelist=req.name.split('.')
      specifier=FS.folder_structure[ltype]
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
for ltype in FS.folder_structure.keys():
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
    FS.folder_structure.update(kwargs)

yaml_classes=[DataFile,UpdateFolderStructure]+fs_locators

