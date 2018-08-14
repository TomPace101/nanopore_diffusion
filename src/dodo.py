"""Doit file for running requests

Basic Usage:
doit
For more information, try
doit help

To specify an alternate control file,
doit control=controlfile"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages
from doit import get_var

#Local
from requesthandler.filepath import Path
import requesthandler.requestfile as requestfile
import requesthandler.locators as locators

#Constants
default_controlfile = locators.datafolder / 'control.yaml'

#Get the controlfile as a Path
controlfile=get_var('control',default_controlfile)
controlfile=Path(controlfile)

#Initialize a RequestFileRequest
req=requestfile.RequestFileRequest(requestfile=controlfile)

def task_all():
  return req.all_tasks()
