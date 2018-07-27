"""Define the base class for all Requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from itertools import chain
import json

#Site packages
try:
  from doit.tools import config_changed
except ImportError:
  #dummy version of config_changed to be used when doit is not available
  def config_changed(arg):
    return arg

class Request(object):
  """Base class for all requests.
  
  Frequently, derived classes will add slots for some of the following attributes, which this class supports:

    - _inputfile_attrs: list of attributes containing input file paths
    - _more_inputfiles: list of additional input file paths not contained in other attributes
    - _outputfile_attrs: list of attributes containing output file paths
    - _more_outputfiles: list of additional output file paths not contained in other attributes
    - _folders: dictionary {attribute name: folder path} used to find the folders for prepending to file names specified in other attributes
    - _required_attrs: list of attribute names that must be defined when the object is first loaded
    - _config_attrs: list of attribute names that contain the "configuration" of the object, to be included when doing a check for changes to configuration
    - _taskname_src_attr: string that is the name of the attribute that contains the task name for the object
    - _child_attrs: list of attributes that contain other Requests
    - _child_seq_attrs: list of attributes that contain sequences of other Requests"""
  __slots__=() #Needed even if empty: without this, a __dict__ object will be created even though subclasses use __slots__
  def __init__(self,**kwargs):
    ##self.__dict__.update(kwd) #Using __slots__ means there is no __dict__
    #Load the attributes specified
    for k,v in kwargs.items():
      setattr(self,k,v)
    #Check for required attributes that are missing
    if hasattr(self,'_required_attrs'):
      missing = [a for a in self._required_attrs if not hasattr(self,a)]
      assert len(missing)==0, "%s missing required attributes: %s"%(type(self),missing)
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
    for attr in self._all_slots:
      itm=getattr(self,attr,None)
      if hasattr(itm,'to_dict') and callable(itm.to_dict):
        d[attr]=itm.to_dict()
      else:
        d[attr]=itm
    return d
  def all_children(self):
    """Generator yielding all the children of this Request
    
    First, yields entries in _child_attrs.
    Then, yields entries in _child_seq_attrs."""
    for cname in getattr(self,'_child_attrs',[]):
      yield getattr(self,cname)
    for lname in gettatr(self,'_child_seq_attrs',[]):
      itr=getattr(self,lname)
      for req in itr:
        yield req
  def all_tasks(self):
    """Generator yielding task definitions from this Request and all child Requests"""
    selftask=self.task_definition
    if selftask is not None:
      yield selftask
    for req in self.all_children():
      for td in req.all_tasks():
        yield td
  @property
  def task_definition(self):
    """A doit task definition dictionary appropriate for this Request
    
    No task is returned if the taskname is None.
    
    To get requests from this task and its children, see yield_tasks."""
    if self.taskname is None:
      return None
    else:
      return {'name': self.taskname,
       'file_dep': self.inputfiles,
       'uptodate': [config_changed(self.config)],
       'targets': self.outputfiles,
       'actions': [(self.run,)]}
  @property
  def config(self):
    """A string representing the configuration of the object, suitable for use by doit.tools.config_changed."""
    # return(str(self.config_dict))
    return(json.dumps(self.config_dict,sort_keys=True))
  @property
  def config_dict(self):
    d=self.to_dict()
    cd=dict([(k,v) for k,v in d.items() if k in self._config_attrs])
    #Don't add configuration of children: look at the output files they generate instead
    #(otherwise changes will cascade even if output files are unaltered)
    ##for childattr in getattr(self,'_child_attrs',[]):
    ##  cd[childattr]=getattr(self,childattr).config
    return cd
  @property
  def taskname(self):
    """A string representing the task name in doit"""
    attrname=getattr(self,_taskname_src_attr,None)
    if attrname is None:
      return None
    else:
      return getattr(self,self._taskname_src_attr)
  @property
  def inputfiles(self):
    """A list of all the inputfiles related to this request"""
    return self._compile_file_list('_inputfile_attrs','_more_inputfiles','inputfiles')
  @property
  def outputfiles(self):
    """A list of all the outputfiles related to this request"""
    return self._compile_file_list('_outputfile_attrs','_more_outputfiles','outputfiles')
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
