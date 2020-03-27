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
MAX_BUFFER_MESSAGES = 100
logging.TIMING = 25

#Classes for logging

class LogRecord(object):
  def __init__(self,name,level,msg,args=None,**kwargs):
    #Get creation time
    self.created = datetime.now()
    #Usual parameters, like the stdlib version
    self.name = name
    self.msg = msg
    self.args = args
    self.levelname = logging.getLevelName(level)
    self.levelno = level
    #Store the additional parameters
    self.parameters=kwargs
  @property
  def message(self):
    msg = str(self.msg)
    if self.args is not None:
      msg = msg % self.args
    return msg

class Logger(logging.Logger):
  def __init__(self, name, level=logging.NOTSET):
    #Initialization from base class
    super(Logger, self).__init__(name,level)
    #Initialize dictionary of timers
    self.timers={}
  # def wrapFindCaller(self,stack_info=False,stacklevel=1,exc_info=None):
  #   #Copied from the python standard library source, in _log
  #   #Won't work yet
  #   sinfo = None
  #   if _srcfile: #TODO: this is set in the stdlib logging module, and used by findCaller, so it probably needs modification somehow
  #     #IronPython doesn't track Python frames, so findCaller raises an
  #     #exception on some versions of IronPython. We trap it here so that
  #     #IronPython can use logging.
  #     try:
  #       fn, lno, func, sinfo = self.findCaller(stack_info, stacklevel)
  #     except ValueError: # pragma: no cover
  #       fn, lno, func = "(unknown file)", 0, "(unknown function)"
  #   else: # pragma: no cover
  #     fn, lno, func = "(unknown file)", 0, "(unknown function)"
  #   if exc_info:
  #     if isinstance(exc_info, BaseException):
  #       exc_info = (type(exc_info), exc_info, exc_info.__traceback__)
  #     elif not isinstance(exc_info, tuple):
  #       exc_info = sys.exc_info()
  #   return fn,lno,exc_info,func,extra,sinfo
  def _log(self, level, msg, args=None, **kwargs): #TODO: arguments would need to be modified to let wrapFindCaller work
    record = logging._logRecordFactory(self.name, level, msg, args, **kwargs)
    self.handle(record)
  def startTimer(self,timername):
    self.timers[timername]=timing.Timer()
    self.log(logging.TIMING,"Starting new timer.",timername=timername)
  def splitTimer(self,timername):
    delta_str=self.timers[timername].split()
    delta_totalsec=self.timers[timername].lap.total_seconds()
    self.log(logging.TIMING,"Reporting elapsed time on timer.",
            timer_name=timername,
            elapsed_sec=delta_totalsec,
            elapsed=delta_str)
  def stopTimer(self,timername):
    delta_str=self.timers[timername].stop()
    delta_totalsec=self.timers[timername].delta.total_seconds()
    self.log(logging.TIMING,"Stopping timer.",
            timer_name=timername,
            elapsed_sec=delta_totalsec,
            elapsed=delta_str)

class RootLogger(Logger):
  def __init__(self,level):
    #Create a logger named "root"
    Logger.__init__(self,"root",level)
    #Set up a BufferHandler for messages while other handlers are being set up
    self.bufferhandler=BufferHandler()
    self.bufferhandler.name="Root Buffer Handler"
    self.addHandler(self.bufferhandler,skipbuffer=True)
  def addHandler(self, hdlr, skipbuffer=False):
    #Use method from base class
    super(Logger, self).addHandler(hdlr)
    if not skipbuffer:
      #Send the buffered messages to the new handler
      self.bufferhandler.dump_to(hdlr)
      #Note how many messsages were sent to the new handler
      self.info("Adding new handler, and sending it previous messages.",
                handler_type=type(hdlr).__name__,
                handler_name=hdlr.name,
                num_messages_output=self.bufferhandler.num_records,
                missing_messages=self.bufferhandler.num_deleted)

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
    od['timestamp']=timing.timestamp(record.created)
    od['area']=record.name
    od['level']=record.levelname
    od['message']=record.message
    if getattr(record,'parameters',None):
      od['parameters']=record.parameters
    return [od]

class BufferHandler(logging.Handler):
  """A handler that just keeps a limited buffer of log records
  
  Note that no formatter is needed, because the records are stored directly."""
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
  def num_records(self):
    return len(self.buffer)
  @property
  def num_deleted(self):
    return self.total_records-len(self.buffer)
  def dump_to(self,other):
    """Ask another handler to handle all the buffered records"""
    for record in self.buffer:
      other.handle(record)

#Introduce the "timing" log level
logging.addLevelName(logging.TIMING,"TIMING")
#Tell the stdlib logging module about these new classes
logging.setLogRecordFactory(LogRecord)
logging.setLoggerClass(Logger)
#Set up the root logger
root=RootLogger(logging.WARNING)
Logger.root=root
#Set up the Manager
Logger.manager = logging.Manager(Logger.root)
Logger.manager.setLoggerClass(Logger)

def getLogger(name=None):
  #Redo the module-level function from stdlib
  if name:
    return Logger.manager.getLogger(name)
  else:
    return root

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
