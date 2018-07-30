"""Define data locators"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#Local
import filepath
import folderstructure as FS

class LocatorBase(object):
  """Dummy class, just to make type checking easier"""
  pass

class InDataFolder(LocatorBase):
  """File location relative to ``datafolder``
  
  The location of ``datafolder`` is specified by the environment variable ``DATALOC``.
  If this environment variable does not exist, a default value is provided.
  See ``folderstructure.py``."""
  parentpath=FS.datafolder
  def __init__(self,*args,**kwargs):
    self.subpath=filepath.Path(*args,**kwargs)
  def path(self,req):
    return self.parentpath / self.subpath

class InRequestFolder(LocatorBase):
  """File in location determined by the Request that needs it"""
  def __init__(self,*args,**kwargs):
    ##TODO
    pass
  def path(self,req):
    ##TODO
    pass

yaml_classes=[InDataFolder]

