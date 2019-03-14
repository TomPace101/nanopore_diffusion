"""Interface for dealing with yaml-loadable classes."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
import io

#Site packages
from ruamel.yaml import YAML
yaml=YAML(typ="safe", pure=True)

#This package

#Dictionary of all registered classes, by their names
all_registered=odict()

#Placeholder for tracking file being loaded
now_loading=[]

def register_classes(class_list):
  """Register the classes that might be loaded from a yaml file, from a list of classes
  
  Arguments:
  
    - class_list = sequence of classes, as classes"""
  for yclass in class_list:
    yaml.register_class(yclass)
    all_registered[yclass.__name__]=yclass

def newloader(yfile=None):
  """For when you need to load more than one yaml file at a time.
  
  Arguments:
  
    - yfile = optional name of file the loader will be used for.
  
  The returned object is an instance of YAML,
  which is not actually called a loader, but I can't figure out the actual name."""
  if yfile is not None:
    global now_loading
    now_loading.append(yfile)
  yy=YAML(typ="safe", pure=True)
  for yclass in all_registered.values():
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

def readfile(fpath,multidoc=False):
  """Read from a file, assuming others may be read at the same time."""
  with open(str(fpath),'r') as fp:
    dat=fp.read()
  yaml=newloader(fpath)
  if multidoc:
    allobj=yaml.load_all(dat)
  else:
    allobj=yaml.load(dat)
  filedone()
  return allobj

def writefile(obj,fpath):
  """Write object to yaml file, overwriting"""
  with open(str(fpath),'w') as fp:
    yaml.dump(obj,fp)
  return
  