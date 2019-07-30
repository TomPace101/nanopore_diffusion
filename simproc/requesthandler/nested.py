"""Get and set values for not just simple attributes or keys, but dotted paths of attributes/keys"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict

#This package
from . import yaml_manager

def to_sequence(dpath):
  if isinstance(dpath,str):
    seq=dpath.split('.')
  else:
    seq=dpath
  return seq

def drill_down(obj,seq):
  head=seq[:-1]
  tail=seq[-1]
  if hasattr(obj,'get_nested'):
    parent=obj.get_nested(head)
  else:
    parent=get_nested(obj,head)
  return parent,tail

def get_nested(obj,dpath):
  """Return the value from the specified attribute/key/index path
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
        The path may contain attributes and dictionary keys, with no need to distinguish between them.
        List indices are also allowed, but only for sequence arguments,
        as strings are not cast to other data types.
  
  Returns the requested data."""
  nxt=obj
  seq=to_sequence(dpath)
  for name in seq:
    if isinstance(name,str) and hasattr(nxt,name):
      nxt = getattr(nxt,name)
    else:
      try:
        #Note that this will work for both lists and dictionaries
        nxt=nxt.__getitem__(name)
      except Exception as the_err:
        raise KeyError('Invalid path %s: No attribute, key, or index %s'%(dpath,name)) from the_err
  return nxt

def get_nested_default(obj,dpath,default=None):
  """Just like get_nested, but return a default object if any part of the path is invalid"""
  try:
    return get_nested(obj,dpath)
  except:
    return default

def set_nested(obj,dpath,val):
  """Set the value at the specified attribute/key/index path
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
    - val = value to assign
  
  No return value."""
  seq=to_sequence(dpath)
  parent,tail=drill_down(obj,seq)
  if hasattr(parent,'__setitem__'):
    parent.__setitem__(tail,val)
  else:
    setattr(parent,tail,val)
  return

def new_odict(obj,dpath):
  """Create a new OrderedDict instance suitable for later use with set_nested
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
  
  No return value."""
  seq=to_sequence(dpath)
  parent,tail=drill_down(obj,seq)
  setattr(parent,tail,odict())
  return

class WithNested(object):
  """A very basic class that loads and sets nested attributes, and supports reading and writing itself to/from yaml"""
  def __init__(self,**kwargs):
    #Load the attributes specified
    for k,v in kwargs.items():
      setattr(self,k,v)
  def to_dict(self):
    """Return a dictionary with all the object's attributes.

    Note that changes to this dictionary will not affect the object.

    No arguments.

    Returns the dictionary."""
    d={}
    for attr,itm in self.__dict__.items():
      if hasattr(itm,'to_dict') and callable(itm.to_dict):
        d[attr]=itm.to_dict()
      else:
        d[attr]=itm
    return d
  def __getstate__(self):
    """Used for pickling, and converting to yaml"""
    return self.to_dict()
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)
  def get_nested(self,dpath):
    """Return the value from the specified attribute/key/index path"""
    return get_nested(self,dpath)
  def get_nested_default(self,dpath,default=None):
    """Return the value from the specified attribute/key/index path, with default value if not present"""
    return get_nested_default(self,dpath,default)
  def set_nested(self,dpath,val):
    """Set the value at the specified attribute/key/index path"""
    set_nested(self,dpath,val)
    return
  def new_odict(self,dpath):
    """Set up a new OrderedDict for later use with set_nested"""
    new_odict(self,dpath)
    return

class Stored(object):
  """Provide the location of a value stored elsewhere in memory
  
    Attributes:
    
      - attrloc = basicaly a dpath argument for use with get_nested and set_nested
  """
  def __init__(self,attrloc):
    self.attrloc=attrloc
  @classmethod
  def from_yaml(cls, constructor, node):
    return cls(node.value)
  @classmethod
  def to_yaml(cls, representer, node):
    return representer.represent_scalar("!"+cls.__name__,str(node.subpath))

#Register for loading from yaml
yaml_manager.register_classes([Stored])
