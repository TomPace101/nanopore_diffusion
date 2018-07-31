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
    """Generic locator class, to be created from the factory function"""
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

yaml_classes=[DataFile]+fs_locators

