"""Doit file for running requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages

#Local
import requestfile

#Constants
controlfile = FS.datafolder / 'control.yaml'

#Initialize a RequestFileRequest
req=requestfile.RequestFileRequest(requestfile=controlfile.fullpath)

def task_all():
  return req.all_tasks()
