"""Get and set values for not just simple attributes or keys, but dotted paths of attributes/keys"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict

def to_sequence(dpath):
  if isinstance(dpath,str):
    seq=dpath.split('.')
  else:
    seq=dpath
  return seq

def drill_down(obj,seq):
  head=seq[:-1]
  tail=seq[-1]
  parent=obj.get_nested(head)
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
      except:
        raise KeyError('Invalid path %s: No attribute, key, or index %s'%(dpath,name))
  return nxt

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
