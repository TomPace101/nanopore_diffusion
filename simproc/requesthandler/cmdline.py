"""Command-line support for running requests, and returning requests as doit tasks"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import doctest
import importlib
import sys

#Site packages
from doit import get_var

#Local
from . import *

#Constants
default_controlfile = locators.DATAFOLDER / 'control.yaml'

def run_validation(verbose=False):
  """Run the validation tests"""
  reslist=[]
  for m in doctest_modules:
    print("Testing module %s."%m.__name__)
    reslist.append(doctest.testmod(m,verbose=verbose))
    print("---")
  for fpath in doctest_files:
    print("Running tests in %s."%fpath)
    reslist.append(doctest.testfile(fpath,module_relative=False,verbose=verbose)) #To make sure the file can be found when the package is zipped, we have already found its absolute path
    print("---")
  fails,atts=[sum(l) for l in zip(*reslist)]
  print("Passed %d/%d total"%(atts-fails,atts))

def run():
  """Run request files from the command line"""
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('requestfile',nargs="*",help="Path to file containing the request(s) to run. Multiple request files may be listed.")
  parser.add_argument('--verbose',action='store_true',help="Provide verbose output where appropriate.")
  parser.add_argument('--tasks_only',action='store_true',help="Run only requests that define tasks.")
  parser.add_argument('--modules',nargs="+",metavar="MODULE",help="Additional python modules defining classes loadable from yaml input")
  parser.add_argument('--validate',action='store_true',help="Perform validation. If requestfiles are also listed, validation is run first.")
  #TODO: allow selecting a subset of the requests?
  cmdline=parser.parse_args()
  
  #run validation if requested
  if cmdline.validate:
    run_validation(verbose=cmdline.verbose)

  #Confirm that specified request file(s) exist(s)
  file_list=[filepath.Path(rf,isFile=True) for rf in cmdline.requestfile]
  nonexist=[rf for rf in file_list if not rf.exists()]
  assert len(nonexist)==0, "Could not find specified request file(s) %s"%str([rf.fullpath for rf in nonexist])

  #Load the requested modules
  if cmdline.modules is not None:
    customization.load_modules(cmdline.modules)

  #Initialize a RequestFileListRequest
  req=requestfile.RequestFileListRequest(requestfiles=file_list)

  #Run
  if cmdline.tasks_only:
    for treq in req.all_task_requests():
      if cmdline.verbose:
        print(treq.name)
      treq.run()
  else:
    req.run()

def yield_doit_tasks():
  """Create task definitions for doit
  
  Usage:
    ```def task_all():
        return yield_doit_tasks()```
  """
  #Get the controlfile as a Path
  controlfile=get_var('control',default_controlfile)
  controlfile=filepath.Path(controlfile)

  #Initialize a RequestFileRequest
  req=requestfile.RequestFileRequest(requestfile=controlfile)

  #Return the task generator
  return req.all_tasks()

