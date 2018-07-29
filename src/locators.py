"""Define data locators"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#Local
import filepath
import folderstructure as FS

class DataFolderFile(object):
  """File location relative to ``datafolder``
  
  The location of ``datafolder`` is specified by the environment variable ``DATALOC``.
  If this environment variable does not exist, a default value is provided.
  See ``folderstructure.py``."""
  parentpath=FS.datafolder
  def __init__(self,*args,**kwargs):
    self.subpath=filepath.Path(*args,**kwargs)
  @property
  def path(self):
    return self.parentpath / self.subpath

class RequestFile(object):
  """File in location determined by the Request that needs it"""
  ##TODO: how does the request provide this information?
  def __init__(self,*args,**kwargs):
    ##TODO
    pass
  @property
  def path(self):
    ##TODO
    pass

yaml_classes=[DataFolderFile]

