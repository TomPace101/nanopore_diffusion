"""Functions and classes relevant for implementing customization"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import sys
import types

#Site packages

#This package
from . import filepath
from . import locators
from . import request

#Locator for module files
locators.folder_structure.update(modulefile=['modules'])

def load_modules(module_name_list):
  """Load the specified list of modules.
  
  Each module must be in a location already reachable, as with an import statement
  
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

def load_module_from_path(fpath):
  """Load the python module at the requested location
  
  The module is NOT added to ``sys.modules``.
  
  Arguments:
  
    - fpath = path to the module file, as a filepath.Path instance
  
  Returns:
  
    - the loaded module"""
  module_name = fpath.stem
  file_path = fpath.fullpath
  #Method for python between 2 and 3.5
  if sys.version_info.major == 2 or sys.version_info.minor < 5:
    import imp
    module = imp.load_source(module_name, file_path)
  #Method for python >= 3.5
  else:
    import importlib.util
    spec = importlib.util.spec_from_file_location(module_name, file_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
  return module

_CustomizableRequest_props_schema_yaml="""#CustomizableRequest
modules:
  anyOf:
    - {type: 'null'}
    - type: array
      items: {type: pathlike}
methods:
  anyOf:
    - {type: 'null'}
    - type: array
      items: {type: string}
    - type: object
      additionalProperties: {type: string}
initializations:
  anyOf:
    - {type: 'null'}
    - {type: object}
extra:
  anyOf:
    - {type: 'null'}
    - {type: object}
_custom_methods:
  anyOf:
    - {type: 'null'}
    - {type: array}"""

class CustomizableRequest(request.Request):
  """A request that can monkey-patch itself
  
  User-defined attributes:
  
    - modules = a sequence of module file paths to be imported.
      If the module contains the variable ``request_methods``,
      as a list of functions (not function names), then those functions will be
      available for binding as methods of the request,
      by using the ``methods`` attribute.
      As such, those functions should have a first argument of ``self``,
      which will be the request to which the method has been bound.

    - methods = dictionary {method_name: function_name} or list (method_name=function_name) of methods to bind.
      The function_name must be for a function listed in the ``request_methods`` variable
      of exactly one of the modules listed in ``modules``.
      This is intended to allow a module to define multiple options for a given method,
      with the request selecting the one it wants.

    - initializations = a dictionary {module name: {variable: value}}

      Upon loading the module, the values in this dictionary will be passed as keyword arguments
        to the function `initialize_module`, if present, within the module.
      Modules listed here but not in `modules` are silently ignored.

    - extra = dictionary {additional request parameters: assigned value}
    
  Calculated attributes:
  
    - _more_inputfiles: the list of module files is added here
    - _custom_methods: list of the names of the custom methods added to the instance"""
  _props_schema=request.make_schema(_CustomizableRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(CustomizableRequest, self).__init__(**kwargs)
    #Read the customization attributes, allowing any to be missing
    self.resolve_locators_in(['modules'],getattr(self,'name',''))
    modules=getattr(self,'modules',[])
    methods=getattr(self,'methods',[])
    initializations=getattr(self,'initializations',{})
    extra=getattr(self,'extra',{})
    #Assign extra attributes
    for k,v in getattr(self,'extra',{}).items():
      setattr(self,k,v)
    #Module files are file dependencies
    self._more_inputfiles=getattr(self,'_more_inputfiles',[])
    self._more_inputfiles+=modules
    #Load modules
    function_name_to_module={} #dictionary mapping addable functions to their home modules
    for modpath in modules:
      #Load module
      modpath=filepath.Path(modpath,isFile=True)
      themod=load_module_from_path(modpath)
      modname = themod.__name__
      #Intialize, if requested
      kwargs = initializations.get(modname,None)
      if (kwargs is not None) and hasattr(themod,'initialize_module'):
        themod.initialize_module(**kwargs)
      #Track names of addable methods in this module
      for function_obj in getattr(themod,'request_methods',[]):
        function_name_to_module[function_obj.__name__]=themod
    #Convert list of methods into dictionary with key=value
    if isinstance(methods,list):
      methods=dict([(k,k) for k in methods])
    #Bind methods
    self._custom_methods=[]
    for method_name, function_name in methods.items():
      assert function_name in function_name_to_module.keys(), "Unable to find function `%s` in any of these modules: %s"%(function_name,modules)
      themod = function_name_to_module[function_name]
      function_obj = getattr(themod,function_name)
      assert isinstance(function_obj,types.FunctionType), "%s in module %s is not a function"%(function_name,themod.__name__)
      self._custom_methods.append(function_name)
      setattr(self,method_name,types.MethodType(function_obj,self)) #Must type cast to MethodType in order to get implicit first argument `self`
  def validate(self):
    d=self.to_dict()
    #Remove attributes added by customization: these aren't validated
    for attr in getattr(self,'_custom_methods',[]):
      o=d.pop(attr)
    for attr in getattr(self,'extra',{}).keys():
      o=d.pop(attr)
    self.validate_kwargs(**d)


#Convenience function for schema updates
make_schema=request.create_schema_updater(CustomizableRequest)
