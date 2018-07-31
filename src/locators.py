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

# class DataFile(filepath.Path):
#   """Data file location relative to ``datafolder``
#   
#   The location of ``datafolder`` is specified by the environment variable ``DATALOC``.
#   If this environment variable does not exist, a default value is provided.
#   See ``folderstructure.py``."""
#   parentpath=FS.datafolder
#   def __new__(cls, *args, **kwargs):
#     subpath=filepath.Path(*args,**kwargs)
#     self = cls.parentpath / subpath
#     return self


yaml_classes=[DataFile]

