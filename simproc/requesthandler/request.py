"""Define the base class for all Requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages
try:
  from doit.tools import config_changed
except ImportError:
  #dummy version of config_changed to be used when doit is not available
  def config_changed(arg):
    return arg

#This package
from . import filepath
from . import yaml_manager
from . import schema
from . import nested

#Validation schema for Request
#Note that validation is not applied to class attributes,
#but some of these could be instance attributes in a subclass
_Request_props_schema_yaml="""#Request
name:
  anyOf:
    - {type: 'null'}
    - {type: string}
_self_task: {type: boolean}
_inputfile_attrs:
  type: array
  items: {type: string}
_more_inputfiles:
  type: array
  items: {type: pathlike}
_outputfile_attrs:
  type: array
  items: {type: string}
_more_outputfiles:
  type: array
  items: {type: pathlike}
_config_attrs:
  type: array
  items: {type: string}
_child_attrs:
  type: array
  items: {type: string}
_child_seq_attrs:
  type: array
  items: {type: string}
_uptodate: {type: boolean}
"""

class Request(schema.SelfValidating):
  """Base class for all requests. Abstract only, not really meant to be instantiated.
  
  All requests have the abilities to:
  
    1) Run themselves. (see the ``run`` method)
    2) Provide ``doit`` task definitions for themselves and any sub-requests. (see the ``all_tasks`` method)
    3) Validate their own configuration. (see the ``validate`` classmethod inherited from ``schema.SelfValidating``)
  
  Note that validation is performed as part of the construction of a request.
  You cannot create incomplete requests that fail validation, and then modify them to be valid.
  
  User-Provided attributes:
  
    - name: a globally unique string identifying the request
    
        Some subclasses may require this to be defined, others may not.
        If defined, the request will be added to an ordered dictionary of requests by name: {request.name: request, ...}
    
  Frequently, derived classes will add member or class attributes for some of the following attributes, which this class supports:

    - _self_task: boolean, True if request itself defines a task, False if only tasks come from children. If False, children may still define tasks. Defaults to False if not defined.
    - _inputfile_attrs: list of attributes containing input file paths
    - _more_inputfiles: list of additional input file paths not contained in other attributes
    - _outputfile_attrs: list of attributes containing output file paths
    - _more_outputfiles: list of additional output file paths not contained in other attributes
    - _required_attrs: list of attribute names that must be defined for the request to be valid
    - _config_attrs: list of attribute names that contain the "configuration" of the object, to be included when doing a check for changes to configuration
    - _child_attrs: list of attributes that contain other Requests
    - _child_seq_attrs: list of attributes that contain sequences of other Requests
    - _props_schema: jsonschema used to validate request configuration, as a dictionary
        The schema is for the 'properties' element only. The rest is provided internally.
    - _uptodate: boolean, True if up-to-date, False otherwise
      Subclasses may replace the attribute with a property,
      and instances may define their own attribute to override as well.

  Subclasses which return their own doit tasks MUST do the following:
  
    - set _self_task to True
    - include 'name' in _required_attrs
    - define _config_attrs
    - have a `run` method which does NOT return anything
    - provide all their input and output files, which may come from locators
        - input files are specified by _inputfiles_attrs and _more_inputfiles
        - output files are specified by _outputfiles_attrs and _more_outputfiles"""
  _props_schema=schema.SelfValidating.update_props_schema(_Request_props_schema_yaml)
  _uptodate=True
  def get_stored(self,candidate):
    """Get the value from the specified storage location"""
    if isinstance(candidate, nested.Stored):
      return self.get_nested(candidate.attrloc)
    else:
      return candidate
  def render(self,loc):
    """Render a locator to a Path instance
    
    If the provided object is not a locator, it is returned unchanged"""
    if hasattr(loc,'path'):
      reqname=getattr(self,'name','')
      fpath=loc.path(reqname)
    else:
      fpath=filepath.Path(loc)
    return fpath
  def renderstr(self,loc):
    """Convenience function for str(self.render(loc))"""
    return str(self.render(loc))
  def pre_run(self):
    """Steps commonly taken just before execution
    
    These steps include final checks and routine preparations
    needed before the actual execution of the request.
    
    This method is meant to be called by run(),
    but this is not required."""
    #Confirm validation
    self.validate()
    #Confirm that required input files exist
    ## self.confirm_inputfiles() #This won't work because not all input files are required.
    #Create directories for output files if necessary
    self.assure_output_dirs()
    return
  def run(self):
    """Run all child requests
    
    Base classes that implement their own tasks should generally override this method."""
    #Final checks and preparatory steps
    self.pre_run()
    #We only need to call run() on the immediate children.
    #Children with their own children will do the same.
    for req in self.all_children():
      req.run()
  @classmethod
  def complete(cls,**kwargs):
    """Convenience function to set up and run the request.

    Arguments:

      - \*\*kwargs to be passed to the the request class __init__"""
    obj=cls(**kwargs)
    obj.run()
    return obj
  def all_children(self):
    """Generator yielding all the immediate children of this Request
    
    First, yields entries in _child_attrs.
    Then, yields entries in _child_seq_attrs.
    
    Yields only immediate children, not grandchildren, etc.:
    see ``recursive_children`` to go deeper."""
    for cname in getattr(self,'_child_attrs',[]):
      yield getattr(self,cname)
    for lname in getattr(self,'_child_seq_attrs',[]):
      itr=getattr(self,lname)
      for req in itr:
        yield req
  def recursive_children(self):
    """Generator recursively yielding all the children of this Request, and their children, etc."""
    for ch in self.all_children():
      yield ch
      for gch in ch.recursive_children():
        yield gch
  @property
  def config_dict(self):
    d=self.to_dict()
    cd=dict([(k,v) for k,v in d.items() if k in self._config_attrs])
    #Don't add configuration of children: look at the output files they generate instead
    #(otherwise changes will cascade even if output files are unaltered)
    return cd
  @property
  def config(self):
    """A string representing the configuration of the object, suitable for use by doit.tools.config_changed."""
    # return(str(self.config_dict))
    return(yaml_manager.writestring(self.config_dict))
  def _compile_file_list(self,attrs_list_attr,files_list_attr,child_attr=None):
    """Construct a list of files, from the following arguments:
    
    - attrs_list_attr: name of attribute containing a list of attribute names, each of which contains one file path
    - files_list_attr: name of attribute containing a list of file paths
    - child_attr: optional, name of attribute of each child containing their corresponding list
    
    If child_attr is not provided, files from children are not included"""
    #Get files from list of attribute names containing files
    attr_list=getattr(self,attrs_list_attr,[])
    fl = [getattr(self,itm) for itm in attr_list if hasattr(self,itm)]
    #Get files from list of additional files
    fl+=getattr(self,files_list_attr,[])
    #Get files from children, if requested
    if child_attr is not None:
      for child in self.all_children():
        fl += getattr(child,child_attr)
    return fl
  @property
  def inputfiles(self):
    """A list of all the inputfiles related to this request"""
    return self._compile_file_list('_inputfile_attrs','_more_inputfiles')
  @property
  def recursive_inputfiles(self):
    """A list of all the inputfiles related to this request and its children"""
    return self._compile_file_list('_inputfile_attrs','_more_inputfiles','recursive_inputfiles')
  @property
  def outputfiles(self):
    """A list of all the outputfiles related to this request"""
    return self._compile_file_list('_outputfile_attrs','_more_outputfiles')
  @property
  def recursive_outputfiles(self):
    """A list of all the outputfiles related to this request and its children"""
    return self._compile_file_list('_outputfile_attrs','_more_outputfiles','recursive_outputfiles')
  # def confirm_inputfiles(self):
  #   """Issue an error for any input files that do not exist
  # 
  #   For this request only; child requests will need to do this themselves."""
  #   missing=[fpath for fpath in self.inputfiles if not fpath.exists()]
  #   errstring="\n".join(["- "+str(fpath) for fpath in missing])
  #   assert len(missing)==0, "Missing required input files:\n%s\n"%errstring
  #   return
  def assure_output_dirs(self):
    """Create directories for output files if necessary
    
    For this request only; child requests will need to do this themselves."""
    for fpath in self.outputfiles:
      self.render(fpath).assure_dir()
    return
  @property
  def task_definition(self):
    """A doit task definition dictionary appropriate for this Request
    
    To get requests from this task or its children, as appropriate, see all_tasks."""
    return {'name': self.name,
     'file_dep': [self.render(p) for p in self.inputfiles],
     'uptodate': [config_changed(self.config),self._uptodate],
     'targets': [self.render(p) for p in self.outputfiles],
     'actions': [(self.run,)]}
  def all_tasks(self):
    """Generator yielding task definitions from this Request or all child Requests"""
    if getattr(self,'_self_task',False):
      #This request defines a task, so return that
      yield self.task_definition
    else:
      #Yield the tasks of the children
      for req in self.all_children():
        for td in req.all_tasks():
          yield td

#Convenience function for schema updates
make_schema=Request.update_props_schema
#Similar functions would be needed for other request base classes

#jsonschema validator setup
# #For jsonschema version 3
# type_checker = ValidatorClass.TYPE_CHECKER
# #type_checker.redefine(#type name as string, #checking function as callable)
# type_checker=type_checker.redefine("path",lambda chkr,inst: isinstance(inst,filepath.Path))
# type_checker=type_checker.redefine("request",lambda chkr,inst: isinstance(inst,Request))
# ValidatorClass = jsonschema.extend(jsonschema.Draft3Validator, type_checker=type_checker)
#For jsonschema 2.6
schema.extra_types_dict['request']=(Request,)

#Register for loading from yaml
yaml_manager.register_classes([Request])
