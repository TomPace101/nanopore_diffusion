"""Get and set values for not just simple attributes or keys, but dotted paths of attributes/keys"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

def get_nested(obj,dpath):
  """Return the value from the specified attribute/key/index path
  
  Arguments:
  
    - dpath = string describing path to the data, using dot separators, or a sequence
        The path may contain attributes and dictionary keys, with no need to distinguish between them.
        List indices are also allowed.
  
  Returns the requested data."""
  nxt=obj
  if isinstance(dpath,str):
    seq=dpath.split('.')
  else:
    seq=dpath
  for name in seq:
    if hasattr(nxt,name):
      nxt = getattr(nxt,name)
    else:
      try:
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
  if isinstance(dpath,str):
    seq=dpath.split('.')
  else:
    seq=dpath
  head=seq[:-1]
  tail=seq[-1]
  parent=obj.get_nested(head)
  if hasattr(parent,'__setitem__'):
    parent.__setitem__(tail,val)
  else:
    setattr(parent,tail,val)
  return
