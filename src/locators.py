"""Define data locators"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#Local
import filepath
import folderstructure as FS


class RelativeLocation(filepath.Path):
  """File location relative to a class-specific location."""
  parentpath=filepath.Path('.')
  def __new__(cls, *args, **kwargs):
    subpath=filepath.Path(*args,**kwargs)
    self = cls.parentpath / subpath
    return self


class DataFile(RelativeLocation):
  """Data file location relative to ``datafolder``
  
  The location of ``datafolder`` is specified by the environment variable ``DATALOC``.
  If this environment variable does not exist, a default value is provided.
  See ``folderstructure.py``."""
  parentpath=FS.datafolder

class RequestFile(filepath.Path):
  """File in location determined by the Request that needs it"""
  ##TODO: how does the request provide this information?
  #Maybe this should be constructed by a method of the request instead
  @classmethod
  def locate(cls,req):
    """Create a new file location based on the request"""
    pass

yaml_classes=[DataFile]

