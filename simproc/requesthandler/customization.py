"""Functions and classes relevant for implementing customization"""

##TODO
##will this work now?

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
      items:
        anyOf:
          - {type: string}
          - {type: path}
initializations:
  anyOf:
    - {type: 'null'}
    - {type: object}
extra:
  anyOf:
    - {type: 'null'}
    - {type: object}"""

class CustomizableRequest(request.Request):
  """A request that can monkey-patch itself
  
  User-defined attributes:
  
    - modules = a sequence of module file paths to be imported.
      If the module contains the variable ``add_methods``,
      as a list of functions, then those functions will be bound
      as methods of the request.

    - initializations = a dictionary {module name: {variable: value}}

      Upon loading the module, the values in this dictionary will be passed as keyword arguments
        to the function `initialize_module`, if present, within the module.
      Modules listed here but not in `modules` are silently ignored.

    - extra = dictionary {additional request parameters: assigned value}
    
  Calculated attributes:
  
    - _more_inputfiles: the list of module files is added here"""
  _props_schema=request.make_schema(_CustomizableRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(CustomizableRequest, self).__init__(**kwargs)
    #Read the customization attributes, allowing any to be missing
    modules=getattr(self,'modules',[])
    initializations=getattr(self,'initializations',{})
    extra=getattr(self,'extra',{})
    #Module files are file dependencies
    self._more_inputfiles=getattr(self,'_more_inputfiles',[])
    self._more_inputfiles+=modules
    #Bind methods
    for modpath in modules:
      #Load module
      modpath=filepath.Path(modpath,isFile=True)
      themod=load_module_from_path(modpath)
      modname = themod.__name__
      #Intialize, if requested
      kwargs = initializations.get(modname,None)
      if (kwargs is not None) and hasattr(themod,'initialize_module'):
        themod.initialize_module(**kwargs)
      #Assign all module functions as methods of this request
      add_methods=getattr(themod,'add_methods',[])
      mod_contents=dict([(f.__name__,f) for f in add_methods])
      for nm, itm in mod_contents.items():
        if isinstance(itm,types.FunctionType):
          setattr(self,nm,types.MethodType(itm,self)) #Must type cast to MethodType in order to get implicit first argument `self`
      #Assign extra attributes
      for k,v in getattr(self,'extra',{}).items():
        setattr(self,k,v)


#Convenience function for schema updates
make_schema=request.create_schema_updater(CustomizableRequest)
