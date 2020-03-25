"""For creating detailed logs

Each log entry has the following:
- timestamp: a timestamp string
- level: a logging level as defined by the python logging module
- message: a text string
- parameters: a dictionary of values to report

Log entries can also request initializing timers and reporting their current elapsed value.

Not implemented:
- entries: a list of sub-entries
"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import collections
import logging
from datetime import datetime

#Site packages
#(implement a separate yaml instance from the one in yaml_manager, to avoid conflicts)
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap
yaml=YAML(typ="rt",pure=True)
yaml.default_flow_style = False

#This package
from . import filepath
from . import timing
from . import yaml_manager
from . import schema

##TODO: we want the logdir to be relative to the data folder, which means we have to wait until it is defined, right?
##TODO: at what point do we create the log directory if it doesn't exist?

#Constants
MAX_BUFFER_MESSAGES=100

#Classes for logging

class YAMLStreamHandler(logging.StreamHandler):
  """For output of logs to a YAML stream"""
  def emit(self,record):
    try:
      msg=self.format(record)
      stream=self.stream
      yaml.dump(msg,stream)
      self.flush()
    except RecursionError:
      raise
    except Exception:
      self.handleError(record)

class YAMLFormatter(object):
  """Intended only for use with the YAMLStreamHandler"""
  def format(self,record):
    od=CommentedMap()
    od['timestamp']=datetime.fromtimestamp(record.created)
    od['area']=record.name
    od['level']=record.levelname
    od['message']=record.getMessage()
    return [od]

class BufferHandler(logging.Handler):
  """A handler that just keeps a limited buffer of log records"""
  def __init__(self,level=logging.NOTSET):
    #Initialization from base class
    super(BufferHandler, self).__init__(level)
    #Initalize the buffer
    self.buffer=collections.deque(maxlen=MAX_BUFFER_MESSAGES)
    #Initialize record count
    self.total_records=0
  def emit(self,record):
    self.buffer.append(record)
    self.total_records+=1
  @property
  def deleted_records(self):
    return self.total_records-len(self.buffer)
  def dump_to(self,other):
    """Ask another handler to handle all the buffered records"""
    for record in self.buffer:
      other.handle(record)

#Functions for configuring

def find_unique_id(stem,logdir,ext,num_digits,sepchar):
  logdir_path=filepath.Path(logdir) #Is that path relative, or absolute?
  assert logdir_path.is_dir(), "logdir must be a directory"
  existing_files=[c.name for c in logdir_path.iterdir()]
  fname_tmpl=stem+sepchar+"{0:0%dd}"%num_digits
  unid=1
  trial_fname=fname_tmpl.format(unid)
  while trial_fname in existing_files:
    unid+=1
    trial_fname=fname_tmpl.format(unid)
  return trial_fname

def configure_logging(stem="simproc",logdir="logs",ext=".log.yaml",num_digits=3,sepchar="."):
  logfile=find_unique_id(stem,logdir,ext,num_digits,sepchar)
  ##TODO
  return None

#Class for configuring logging from a yaml file

#Validation schema for ConfigLogging
_ConfigLogging_props_schema_yaml="""#ConfigLogging
level: {type: string}
destination: {type: pathlike}
"""

class ConfigLogging(schema.SelfValidating):
  """Class to configure logging parameters from within a yaml file

  This isn't really a class. It's a function.
  But when loading from yaml, there isn't a way to call a function directly.
  You can only load classes.
  Hence, this is an object that simply calls another function when initialized,
  and then does nothing else ever.
  
  This is also not a request: its action is taken at initialization,
  not when requests are run.

  Initialization arguments:
    - level: minimum level for events to log
    - destination: path to the log output file
  """
  _props_schema=schema.SelfValidating.update_props_schema(_ConfigLogging_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RequestFileRequest, self).__init__(**kwargs)
    ##TODO
 
#Register for loading from yaml
yaml_manager.register_classes([ConfigLogging])
