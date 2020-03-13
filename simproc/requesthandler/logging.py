"""For creating detailed logs

Each log entry has the following:
- timestamp: a timestamp string
- level: a logging level as defined by the python logging module
- message: a text string
- entries: a list of sub-entries
"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import logging

#Site packages
#(implement a separate yaml instance from the one in yaml_manager, to avoid conflicts)
from ruamel.yaml import YAML
yaml=YAML(typ="safe", pure=True)

#This package
from . import yaml_manager
from . import schema

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
