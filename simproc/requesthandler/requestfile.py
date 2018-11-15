"For reading requests from yaml files"

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#This package
from . import filepath
from . import yaml_manager
from . import locators
from . import request

_RequestFileRequest_props_schema_yaml="""#RequestFileRequest
requestfile:
  anyOf:
    - {type: string}
    - {type: path}
_children: {type: array}"""

class RequestFileRequest(request.Request):
  """Request to run all the requests in a given file
  
  User-Provided Attributes:
  
    - requestfile: path to the file containing the requests, as instance of filepath.Path or Locator
    
  Calculated Attributes:
  
    - _children: A list storing all child requests"""
  _props_schema=request.make_schema(_RequestFileRequest_props_schema_yaml)
  _required_attrs=['requestfile']
  _child_seq_attrs=['_children']
  ##_inputfile_attrs=['requestfile'] #I can't decide if this should be here or not
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileRequest, self).__init__(**kwargs)
    #Read the file
    rfpath = self.requestfile.fullpath if hasattr(self.requestfile,'fullpath') else self.requestfile
    with open(rfpath,'r') as fp:
      dat=fp.read()
    #Load all objects from yaml
    yaml=yaml_manager.newloader(rfpath)
    allobj=yaml.load_all(dat)
    #Store child objects that are Request subclasses
    self._children=[ch for ch in allobj if isinstance(ch,request.Request)]
    #Loading of yaml file is complete
    yaml_manager.filedone() #Let the manager know we're no longer loading this file

_RequestFileListRequest_props_schema_yaml="""#RequestFileListRequest
requestfiles:
  type: array
  items:
    anyOf:
      - {type: string}
      - {type: path}
_children: {type: array}"""

class RequestFileListRequest(request.Request):
  """Request to run all of the requests in a given list of files
  
  User-Provided Attributes:
  
    - requestfiles: sequence of paths to the request files, each an instance of filepath.Path or Locator
  
  Calculated Attributes:
  
    - _children: A list storing all child requests
    """
  _props_schema=request.make_schema(_RequestFileListRequest_props_schema_yaml)
  _required_attrs=['requestfiles']
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileListRequest, self).__init__(**kwargs)
    #List inputfiles
    ##self._more_inputfiles=self.requestfiles #I can't decide if this should be here or not
    #Each listed file is a RequestFileRequest
    self._children=[RequestFileRequest(requestfile=ch) for ch in self.requestfiles]

#Register locators and default folder structure
locators.folder_structure.update(RequestFile=['requests',0,1,2,3])

#Register for loading from yaml
yaml_manager.register_classes([RequestFileRequest, RequestFileListRequest])
