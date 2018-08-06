"""Interface for dealing with yaml-loadable classes."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages
from ruamel.yaml import YAML
yaml=YAML(typ="safe", pure=True)

#This package

def register_classes(class_list):
  """Register the classes that might be loaded from a yaml file, from a list of classes
  
  Arguments:
  
    - class_list = sequence of classes, as classes"""
  for yclass in class_list:
    yaml.register_class(yclass)

def read(s):
  """Syntactic sugar for yaml.load()"""
  return yaml.load(s)
