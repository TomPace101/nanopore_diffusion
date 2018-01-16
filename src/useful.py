
#Standard library
import argparse
from functools import reduce
from itertools import chain
import operator
import os.path as osp
import pickle

#Site packages
import yaml
try:
  from doit.tools import config_changed
except ImportError:
  #dummy version of config_changed to be used when doit is not available
  def config_changed(arg):
    return arg

#Constants
pickle_protocol = 4 #The newest protocol, requires python 3.4 or above.

#-------------------------------------------------------------------------------
#Get a value from a nested dictionary

def nested_location(obj,loc):
  """For nested dictionary obj, get item at location specified by the sequence loc

  >>> obj={'a':{'b':{'c':11,'d':99}}}
  >>> loc=['a','b','c']
  >>> nested_location(obj,loc)
  11
  """
  loclist=list(loc)
  return reduce(operator.getitem, [obj]+loclist)

def nested_to_str(obj):
  """For nested dictionary obj, return a consistent string"""
  #Get pairs in correct order, with keys converted to strings
  pairs1=[(repr(k),v) for k,v in obj.items()]
  pairs1.sort(key=lambda t: t[0])
  #Get list of sorted key strings alone
  skeys=[p[0] for p in pairs1]
  #Convert values to strings while maintaining order
  svals=[]
  for k,v in pairs1:
    if hasattr(v,'items'):
      svals.append(nested_to_str(v))
    else:
      svals.append(repr(v))
  items_list=[k+': '+v for k,v in zip(skeys,svals)]
  return '{%s}'%', '.join(items_list)

#-------------------------------------------------------------------------------
#Standard work with yaml files

def readyaml(fpath):
  "read object from yaml file"
  with open(fpath,'r') as fp:
    dat=fp.read()
    obj=yaml.load(dat)
  return obj

def readyaml_multidoc(fpath):
  "Return multiple documents from a yaml file as a list"
  with open(fpath,'r') as fp:
    dat=fp.read()
    gen=yaml.load_all(dat)
  return [obj for obj in gen]

def readyaml_gen(fpath):
  "return a generator for the contents of the yaml file"
  with open(fpath,'r') as fp:
    dat=fp.read()
  gen=yaml.load_all(dat)
  return gen

def writeyaml(obj,fpath):
  "write object to yaml file, overwriting"
  with open(fpath,'w') as fp:
    yaml.dump(obj,fp)
  return

#-------------------------------------------------------------------------------
#Standard work with pickle files

def readpickle(fpath):
  "read object from pickle file"
  with open(fpath, 'rb') as fp:
    obj=pickle.load(fp)
  return obj

def writepickle(obj,fpath):
  "write object to pickle file, overwriting"
  with open(fpath, 'wb') as fp:
    pickle.dump(obj,fp,pickle_protocol)
  return

#-------------------------------------------------------------------------------
#ParameterSet

def getbasename(fpath):
  "Return the stem name of the specified file (no directory, no extension)"
  return osp.splitext(osp.basename(fpath))[0]

class ParameterSet:
  """A base class for defining sets of related parameters.
  Each subclass should use __slots__ to define its parameters.
  This is not intended as a method for storing complicated objects;
  all attributes should have values that are numbers, strings, sequences, or dictionaries
  whose items follow the same rules."""
  __slots__=('basename',) #Needed even if empty: without this, a __dict__ object will be created even though subclasses use __slots__
  def __init__(self,**kwd):
    ##self.__dict__.update(kwd) #Using __slots__ means there is no __dict__
    for k,v in kwd.items():
      setattr(self,k,v)
    #Check for required attributes that are missing
    if hasattr(self,'_required_attr'):
      missing = [a for a in self._required_attr if not hasattr(self,a)]
      assert len(missing)==0, "%s missing required attributes: %s"%(type(self),missing)
  @classmethod
  def from_yaml(cls,fpath):
    """Read ParameterSet from a yaml file.
    Arguments:
      fpath = path to the yaml file to read in
        This must be a single-document yaml file,
        and it is assumed to be structured as a single dictionary at the top level.
    Returns:
      pset = a ParameterSet object as defined by the contents of the yaml file"""
    d=readyaml(fpath)
    obj=cls(**d)
    obj.basename=getbasename(fpath)
    return obj
  @classmethod
  def all_from_yaml(cls,fpath):
    """Generator to read a list of ParameterSet objects from a yaml file.
    Arguments:
      fpath = path to the yaml file to read in
        This must be a multi-document yaml file,
        and each document is assumed to be structured as a single dictionary at the top level.
    Each call yields:
      pset = a ParameterSet object as defined by one document in the yaml file"""
    with open(fpath,'r') as fp:
      dat=fp.read()
      gen=yaml.load_all(dat)
    for d in gen:
      obj=cls(**d)
      obj.basename=getbasename(fpath)
      yield obj
  def to_yaml(self,fpath):
    """Write ParameterSet to a yaml file.
    Arguments:
      fpath = path to the yaml file to write
        This will be a single-document yaml file, containing a dictionary (potentially of other dictionaries).
        The file will be overwritten if it already exists.
    No return value."""
    #TODO: consider writing docstring as comments in yaml file
    writeyaml(self.to_dict(),fpath)
    return
  @classmethod
  def from_pickle(cls,fpath):
    """Read ParameterSet from a pickle file.
    Arguments:
      fpath = path to the pickle file to read in
        This file should be a mapping type at the top level.
    Returns:
      pset = a ParameterSet object as defined by the contents of the pickle file."""
    d=readpickle(fpath)
    return cls(**d)
  def to_pickle(self,fpath):
    """Write ParameterSet to a pickle file.
    Arguments:
      fpath = path to the pickle file to write
        The file will be overwritten if it already exists.
    No return value."""
    writepickle(self.to_dict(),fpath)
    return
  @classmethod
  def from_dict(cls,d):
    return cls(**d)
  def _all_slots_iter(self):
    """Return an iterator over all the available slots"""
    return chain.from_iterable(getattr(cls, '__slots__', []) for cls in type(self).__mro__)
  def to_dict(self):
    """Return a dictionary with all the object's attributes.
    Note that changes to this dictionary will not affect the object.
    No arguments.
    Returns the dictionary."""
    return dict([(k,getattr(self,k,None)) for k in self._all_slots_iter()])
  def __repr__(self):
    return nested_to_str(self.to_dict())
  @classmethod
  def from_Namespace(cls,ns):
    return cls(**ns.__dict__)
  def to_Namespace(self):
    """Return an argparse.Namespace object with all the object's attributes.
    Note that changes to the Namespace will not affect the object.
    No arguments.
    Returns the Namespace."""
    return argparse.Namespace(**self.to_dict())
  ##TODO: read and write from ini file
  #Methods supporting task definition for doit
  @property
  def task_definition(self):
    """A doit task definition dictionary appropriate for processing this object."""
    return {'name': self.taskname,
     'file_dep': list(self.inputfiles.values()),
     'uptdodate': config_changed(self.config()),
     'targets': list(self.outputfiles.values()),
     'actions': [(self.run,)]}
  @property
  def config(self):
    """A string representing the configuration of the object, suitable for use by doit.tools.config_changed."""
    d=self.to_dict()
    cd=dict([(k,v) for k,v in d.items() if k in self._config_attr])
    return(str(cd))
  @property
  def taskname(self):
    """A string representing the task name in doit"""
    return getattr(self,self._taskname_src)

#-------------------------------------------------------------------------------
#Common command-line usage

def run_cmd_line(program_description,input_file_description,objtype,process_method="run",other_selection=None,f_args=None):
  """Perform common command line processing.
  Many of the modules use a similar process for running from the command line:
  - read an multi-document input yaml file
  - load each document into an object of some type
  - select specific objects using the "--select" command line argument, and/or other means
  - run a function on each selected object
  This function implements that common process.
  Arguments:
    program_description = string containing the help for the program to run
    input_file_desscription = string containing help for the input file to process
    objtype = type all documents in the input file should be loaded into
      This is assumed to be a subclass of ParameterSet.
      (At a minimum, it must have an all_from_yaml method)
    process_method = optional string, name of object's method to be called
    other_selection = optional dictionary with additional requirements for objects to process
      {attribute name: sequence of allowed values}
    f_args = optional list of positional arguments for calling the process_method
  No return value."""

  #Parse arguments
  parser = argparse.ArgumentParser(description=program_description)
  parser.add_argument('params_yaml', help=input_file_description)
  parser.add_argument('--select',nargs="+",help="""Only process selected elements of the input yaml file.
    This option must be followed by an attribute name, and then a sequence of values for that attribute.
    Only those entries in the yaml file where the attribute matches one of these values will be processed.""")
  cmdline=parser.parse_args()
  #Check that input file exists
  assert osp.isfile(cmdline.params_yaml), "Parameter definition file does not exist: %s"%cmdline.params_yaml
  #Set up dictionary for object selection by attribute
  if other_selection is None:
    requirements = {}
  else:
    requirements=dict(other_selection)
  if cmdline.select is not None:
    requirements[cmdline.select[0]]=cmdline.select[1:]
  #Any other positional arguments?
  if f_args is None:
    f_args = []

  #Read in the yaml file
  allobjs=objtype.all_from_yaml(cmdline.params_yaml)

  #Process documents
  for obj in allobjs:
    #Is this document selected? (if no selection list provided, process all documents)
    selected=True
    for selattr,allowed in requirements.items():
       selected = selected and (getattr(obj,selattr,None) in allowed)
       if not selected:
         break
    if selected:
      getattr(obj,process_method)(*f_args)
