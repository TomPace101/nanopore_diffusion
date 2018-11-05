"""Interface for dealing with yaml-loadable classes."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import io

#Site packages
from ruamel.yaml import YAML
yaml=YAML(typ="safe", pure=True)

#This package

#List of all registered class names
all_registered=[]
all_registered_classes=[]

#Placeholder for tracking file being loaded
now_loading=[]

def register_classes(class_list):
  """Register the classes that might be loaded from a yaml file, from a list of classes
  
  Arguments:
  
    - class_list = sequence of classes, as classes"""
  for yclass in class_list:
    yaml.register_class(yclass)
    all_registered.append(yclass.__name__)
    all_registered_classes.append(yclass)

def newloader(yfile=None):
  """For when you need to load more than one yaml file at a time.
  
  Arguments:
  
    - yfile = optional name of file the loader will be used for.
  
  The returned object is an instance of YAML,
  which is not actually called a loader, but I can't the actual name."""
  if yfile is not None:
    global now_loading
    now_loading.append(yfile)
  yy=YAML(typ="safe", pure=True)
  for yclass in all_registered_classes:
    yy.register_class(yclass)
  return yy

def filedone():
  """Call to indicate that loading of a file is complete"""
  global now_loading
  now_loading.pop()

def read(s):
  """Syntactic sugar for yaml.load()"""
  return yaml.load(s)

def yamlstring(obj):
  """like json.dumps but for yaml"""
  with io.StringIO() as strm:
    yaml.dump(obj,strm)
    dat=strm.getvalue()
  return dat
