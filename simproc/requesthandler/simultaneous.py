"For running requests simultaneously"

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call
import time

#This package
from . import request
from . import yaml_manager

_SimultaneousRequestQueue_props_schema_yaml="""#SimultaneousRequestQueue
name: {type: string}
num_workers:
  type: integer
  minimum: 1
queue: {type: array}
delay:
  type: number
  minimum: 0
"""


class SimultaneousRequestQueue(request.Request):
  """Run a queue of requests, with more than one allowed to run simultaneously
  
  User-defined attributes:
  
    - num_workers: number of requests to run in parallel
    - queue: sequence of the requests to run
    - delay: optional, time (in seconds), to wait between checks for request completion, defaults to 30"""
  _self_task=False
  _required_attrs=['num_workers','queue']
  _props_schema=request.make_schema(_SimultaneousRequestQueue_props_schema_yaml)
  _child_seq_attrs=['queue']
  def run(self):
    ##TODO
    pass

#Register for loading from yaml
yaml_manager.register_classes([SimultaneousRequestQueue])

