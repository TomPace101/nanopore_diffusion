
#Standard library
import pickle

#Site packages
import yaml

#Constants
pickle_protocol = 4 #The newest protocol, requires python 3.4 or above.

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

class ParameterSet:
  """A base class for defining sets of related parameters.
  Each subclass should use __slots__ to define its parameters.
  This is not intended as a method for storing complicated objects;
  all attributes should have values that are numbers, strings, sequences, or dictionaries
  whose items follow the same rules."""
  __slots__=tuple() #Without this, a __dict__ object will be created even though subclasses use __slots__
  def __init__(self,**kwd):
    ##self.__dict__.update(kwd) #This doesn't seem work well with __slots__
    for k,v in kwd.items():
      setattr(self,k,v)
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
    return cls(**d)
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
    for obj in gen:
      yield cls(**obj)
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
  def to_dict(self):
    """Return a dictionary with all the object's attributes.
    Note that changes to this dictionary will not affect the object.
    No arguments.
    Returns the dictionary."""
    return dict([(k,getattr(self,k,None)) for k in self.__slots__])
  ##TODO: read and write from ini file
