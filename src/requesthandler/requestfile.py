"For reading requests from yaml files"

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#This package
from . import filepath
from . import yaml_manager
from . import request

#object to load/dump yaml
yaml=yaml_manager.yaml

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
  _props_schema=yaml.load(_RequestFileRequest_props_schema_yaml)
  _required_attrs=['requestfile']
  _child_seq_attrs=['_children']
  _inputfile_attrs=['requestfile']
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileRequest, self).__init__(**kwargs)
    #Read the file
    rfpath = self.requestfile.fullpath if hasattr(self.requestfile,'fullpath') else self.requestfile
    with open(rfpath,'r') as fp:
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
  _props_schema=yaml.load(_RequestFileListRequest_props_schema_yaml)
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

#Register for loading from yaml
yaml_manager.register_classes([RequestFileRequest, RequestFileListRequest])
