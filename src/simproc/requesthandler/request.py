"""Define the base class for all Requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from collections import OrderedDict as odict
from itertools import chain

#Site packages
import jsonschema
try:
  from doit.tools import config_changed
except ImportError:
  #dummy version of config_changed to be used when doit is not available
  def config_changed(arg):
    return arg

#This package
from . import filepath
from .yaml_manager import yaml, yamlstring

#Validation partial setup (some setup must wait for Request class to be defined)
ValidatorClass = jsonschema.Draft4Validator
#jsonschema 2.6
extra_types_dict={'path':filepath.Path}

#Dictionary of all loaded requests
all_requests=odict()

def validation_error_string(err):
  "Return a string explaining the given validation error"
  #Get the basic message
  s=err.message
  #Provide a path if necessary
  if len(err.path)>0:
    s="%s: %s"%('.'.join([str(itm) for itm in err.path]),s)
  return s

class Request(object):
  """Base class for all requests. Abstract only, not really meant to be instantiated.
  
  All requests should have the following abilities:
  
    1) Run themselves. (see the ``run`` method)
    2) Provide ``doit`` task definitions for themselves and any sub-requests. (see the ``all_tasks`` method)
    3) Validate their own configuration. (see the ``validate`` classmethod)
  
  Note that validation is performed as part of the construction of a request.
  You cannot create incomplete requests that fail validation, and then modify them to be valid.
  
  User-Provided attributes:
  
    - name: a globally unique string identifying the request
    
        Some subclasses may require this to be defined, others may not.
        If defined, the request will be added to an ordered dictionary of requests by name: {request.name: request, ...}
    
  Frequently, derived classes will add slots or class attributes for some of the following attributes, which this class supports:

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

  Subclasses which return their own doit tasks MUST do the following:
  
    - set _self_task to True
    - include 'name' in _required_attrs
    - define _config_attrs
    - have a `run` method which does NOT return anything
    - provide all their input and output files, which may come from locators
        - input files are specified by _inputfiles_attrs and _more_inputfiles
        - output files are specified by _outputfiles_attrs and _more_outputfiles"""
  __slots__=('name',) #Needed even if empty: without this, a __dict__ object will be created even though subclasses use __slots__
  def __init__(self,**kwargs):
    #Process locators
    for k,v in kwargs.items():
      #If field is a locator, get the Path it returns
      ##TODO: this doesn't catch entries that are lists (check each element)
      ##or nested dictionaries, or lists of dictionaries of lists of ...
      ##how do we make this recursive?
      if hasattr(v,'path'):
        kwargs[k]=v.path(self)
    #Validate kwargs
    if hasattr(self,'_props_schema'):
      self.validate_kwargs(**kwargs)
    #Load the attributes specified
    for k,v in kwargs.items():
      setattr(self,k,v)
    #If named, add to dictionary of named requests
    if hasattr(self,'name'):
      all_requests[self.name]=self
  @classmethod
  def _class_schema(cls):
    """Return the jsonschema validation schema for instances of this class"""
    return {'type':'object',
            'properties':cls._props_schema,
            'required':getattr(cls,'_required_attrs',[]),
            'additionalProperties':False}
  @classmethod
  def validate_kwargs(cls,**kwargs):
    if hasattr(cls,'_props_schema'):
      schema=cls._class_schema()
      validator=ValidatorClass(schema,types=extra_types_dict)
      errlist=["  - %s"%validation_error_string(err) for err in validator.iter_errors(kwargs)]
      if len(errlist)>0:
        #Found errors: raise exception listing them all
        errlist.sort()
        errstr="Errors found in %s.\n"%cls.__name__
        keylist=list(kwargs.keys())
        keylist.sort()
        errstr+='Received arguments:\n'
        errstr+='\n'.join(['  - %s: %s'%(k,kwargs[k]) for k in keylist])
        errstr+='\nErrors:\n'
        errstr+='\n'.join(errlist)
        raise Exception(errstr)
  def validate(self):
    d=self.to_dict()
    self.validate_kwargs(**d)
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)
  def run(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'run' method."%str(type(self)))
  def _all_slots(self):
    """Return an iterator over all the available slots"""
    return chain.from_iterable(getattr(cls, '__slots__', []) for cls in type(self).__mro__)
  def to_dict(self):
    """Return a dictionary with all the object's attributes.

    Note that changes to this dictionary will not affect the object.

    No arguments.

    Returns the dictionary."""
    d={}
    for attr in self._all_slots():
      itm=getattr(self,attr,None)
      if hasattr(itm,'to_dict') and callable(itm.to_dict):
        d[attr]=itm.to_dict()
      else:
        d[attr]=itm
    return d
  def __getstate__(self):
    """Used for pickling, and possibly for converting to yaml"""
    return self.to_dict()
  def all_children(self):
    """Generator yielding all the children of this Request
    
    First, yields entries in _child_attrs.
    Then, yields entries in _child_seq_attrs."""
    for cname in getattr(self,'_child_attrs',[]):
      yield getattr(self,cname)
    for lname in getattr(self,'_child_seq_attrs',[]):
      itr=getattr(self,lname)
      for req in itr:
        yield req
  @property
  def config_dict(self):
    d=self.to_dict()
    cd=dict([(k,v) for k,v in d.items() if k in self._config_attrs])
    #Don't add configuration of children: look at the output files they generate instead
    #(otherwise changes will cascade even if output files are unaltered)
    #All Paths must be converted to strings
    #TODO: it might be preferable to allow Paths to convert themselves to yaml. I couldn't get that to work before, though.
    for k,v in cd.items():
      if isinstance(v,filepath.Path):
        cd[k]=v.fullpath  ##TODO: this means moving/renaming the data folder will show up as a changed configuration for everything
    return cd
  @property
  def config(self):
    """A string representing the configuration of the object, suitable for use by doit.tools.config_changed."""
    # return(str(self.config_dict))
    return(yamlstring(self.config_dict))
  def _compile_file_list(self,attrs_list_attr,files_list_attr,child_attr):
    """Construct a list of files, from the following arguments:
    
    - attrs_list_attr: name of attribute containing a list of attribute names, each of which contains one file path
    - files_list_attr: name of attribute containing a list of file paths
    - child_attr: name of attribute of each child containing their corresponding list"""
    #Get files from list of attribute names containing files
    attr_list=getattr(self,attrs_list_attr,[])
    fl = [itm.fullpath for itm in attr_list]
    #Get files from list of additional files
    fl+=getattr(self,files_list_attr,[])
    #Get files from children
    for child in self.all_children():
      if child.taskname is None: #Children that define their own tasks are responsible for their own file dependencies
        fl += getattr(child,child_attr)
    return fl
  @property
  def inputfiles(self):
    """A list of all the inputfiles related to this request"""
    return self._compile_file_list('_inputfile_attrs','_more_inputfiles','inputfiles')
  @property
  def outputfiles(self):
    """A list of all the outputfiles related to this request"""
    return self._compile_file_list('_outputfile_attrs','_more_outputfiles','outputfiles')
  @property
  def task_definition(self):
    """A doit task definition dictionary appropriate for this Request
    
    No task is returned if the taskname is None.
    
    To get requests from this task and its children, see yield_tasks."""
    return {'name': self.name,
     'file_dep': self.inputfiles,
     'uptodate': [config_changed(self.config)],
     'targets': self.outputfiles,
     'actions': [(self.run,)]}
  def all_tasks(self):
    """Generator yielding task definitions from this Request and all child Requests"""
    if getattr(self,'_self_task',False):
      yield self.task_definition
    for req in self.all_children():
      for td in req.all_tasks():
        yield td

#jsonschema validator setup
# #For jsonschema version 3
# type_checker = ValidatorClass.TYPE_CHECKER
# #type_checker.redefine(#type name as string, #checking function as callable)
# type_checker=type_checker.redefine("path",lambda chkr,inst: isinstance(inst,filepath.Path))
# type_checker=type_checker.redefine("request",lambda chkr,inst: isinstance(inst,Request))
# ValidatorClass = jsonschema.extend(jsonschema.Draft3Validator, type_checker=type_checker)
#For jsonschema 2.6
extra_types_dict['request']=(Request,)

