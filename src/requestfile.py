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
yaml_module_list=['locators','request']

def _direct_register_classes(class_list):
  """Register the classes that might be loaded from a yaml file, from a list of classes
  
  Arguments:
  
    - class_list = sequence of classes, as classes"""
  for yclass in class_list:
    yaml.register_class(yclass)
    # all_classes[yclass.__name__]=yclass

def register_classes(module_list):
  """Register the classes that might be loaded from a yaml file, from a list of module names
  
  Arguments:
  
    - module_list = list of module names, as strings
  
  Each module MUST define the variable ``yaml_classes``,
  as a sequence of classes that can be loaded from yaml file
  
  No return value."""
  for modname in yaml_module_list:
    if modname in sys.modules.keys():
      loaded_module=sys.modules[modname]
    else:
      loaded_module=importlib.import_module(modname)
    class_list=getattr(loaded_module,'yaml_classes',[])
    _direct_register_classes(class_list)

_RequestFileRequest_props_schema_yaml="""#RequestFileRequest
name: {type: string}
requestfile:
  anyOf:
    - {type: string}
    - {type: path}"""

class RequestFileRequest(request.Request):
  """Request to run all the requests in a given file
  
  User-Provided Attributes:
  
    - requestfile: path to the file containing the requests, as instance of filepath.Path or Locator
    
  Calculated Attributes:
  
    - _children: A list storing all child requests"""
  __slots__=('requestfile','_children')
  _props_schema_yaml=_RequestFileRequest_props_schema_yaml
  _required_attrs=['requestfile']
  _child_seq_attrs=['_children']
  _inputfile_attrs=['requestfile']
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileRequest, self).__init__(**kwargs)
    #Read the file
    with open(self.requestfile.fullpath,'r') as fp:
      dat=fp.read()
    #Load all objects from yaml
    allobj=yaml.load_all(dat)
    self._children=[ch for ch in allobj]
  def run(self):
    "Run all the requests listed in the file"
    for req in self.all_children():
      req.run()

_RequestFileListRequest_props_schema_yaml="""#RequestFileListRequest
name: {type: string}
requestfiles:
  type: array
  items:
    anyOf:
      - {type: string}
      - {type: path}"""

class RequestFileListRequest(request.Request):
  """Request to run all of the requests in a given list of files
  
  User-Provided Attributes:
  
    - requestfiles: sequence of paths to the request files, each an instance of filepath.Path or Locator
  
  Calculated Attributes:
  
    - _children: A list storing all child requests
    """
  __slots__=('requestfiles','_children')
  _props_schema_yaml=_RequestFileListRequest_props_schema_yaml
  _required_attrs=['requestfiles']
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
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
