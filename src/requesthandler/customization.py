"""Functions and classes relevant for implementing customization"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import sys

def load_modules(module_name_list):
  """Register the classes that might be loaded from a yaml file, from a list of module names
  
  Arguments:
  
    - module_name_list = list of module names, as strings
  
  Returns:
  
    - loaded_module_list = list of loaded modules, as modules"""
  loaded_module_list=[]
  for modname in module_name_list:
    if modname in sys.modules.keys():
      loaded_module=sys.modules[modname]
    else:
      loaded_module=importlib.import_module(modname)
    loaded_module_list.append(loaded_module)
  return loaded_module_list
