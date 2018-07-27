"For reading requests from yaml files"

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import sys

#Site packages
from ruamel.yaml import YAML
yaml=YAML(typ="safe", pure=True)

#Local
import folderstructure as FS
import filepath
import request

#Complete list of all modules defining classes we want to load from yaml
yaml_module_list=['locators']

def _direct_register_classes(class_list):
  """Register the classes that might be loaded from a yaml file, from a list of classes
  
  Arguments:
  
    - class_list = sequence of classes, as classes"""
  for yclass in class_list:
    yaml.register_class(yclass)

def register_classes(module_list):
  """Register the classes that might be loaded from a yaml file, from a list of module names
  
  Arguments:
  
    - module_list = list of module names, as strings
  
  Each module MUST define the variable ``yaml_classes``,
  as a sequence of classes that can be loaded from yaml file
  
  No return value."""
  for modname in yaml_module_list:
    loaded_module=importlib.import_module(modname)
    class_list=getattr(loaded_module,'yaml_classes',[])
    _direct_register_classes(class_list)

class RequestFileRequest(request.Request):
  """Request to run all the requests in a given file
  
  User-Provided Attributes:
  
    - requestfile: File Locator to the file containing the requests
    
  Calculated Attributes:
  
    - _children: A list storing all child requests"""
  __slots__=('requestfile','_children')
  _required_attrs=['requestfile']
  _child_seq_attrs=['_children']
  _taskname_src_attr=None #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileRequest, self).__init__(**kwargs)
    #Read the file
    with open(self.requestfile.fullpath,'r') as fp:
      dat=fp.read()
    #Load all objects from yaml
    allobj=yaml.load_all(dat)
    #Some of the objects may not be requests; skip those
    self._children=[ch for ch in allobj if isinstance(ch,request.Request)]
  def run(self):
    "Run all the requests listed in the file"
    for req in self.all_children():
      req.run()

class RequestFileListRequest(request.Request):
  """Request to run all of the requests in a given list of files
  
  User-Provided Attributes:
  
    - requestfiles: sequence of paths to the request files
  
  Calculated Attributes:
  
    - _children: A list storing all child requests
    """
  __slots__=('requestfiles','_children')
  _required_attrs=['requestfiles']
  _child_seq_attrs=['_children']
  _taskname_src_attr=None #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileListRequest, self).__init__(**kwargs)
    #Each listed file is a RequestFileRequest
    self._children=[RequestFileRequest(requestfile=ch) for ch in self.requestfiles]
  def run(self):
    "Run all Requests in each listed request file"
    for ch in self.all_children():
      ch.run()

register_classes(yaml_module_list)
yaml_classes=[RequestFileRequest, RequestFileListRequest]
_direct_register_classes(yaml_classes)
