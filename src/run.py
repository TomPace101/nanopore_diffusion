"""Command-line support for running requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import importlib
import sys

#Site packages

#Local
import requestfile

#Handle command-line execution
if __name__ == '__main__':
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('requestfile', help="Path to file containing the request(s) to run")
  parser.add_argument('--modules',nargs="+",help="Additional python modules defining classes loadable from yaml input")
  #TODO: allow selecting a subset of the requests?
  cmdline=parser.parse_args()
  requestfile=filepath.Path(cmdline.requestfile,isFile=True)
  
  #Confirm that specified request file exists
  assert requestfile.exists(), "Could not find specified request file %s"%requestfile.fullpath

  #Add the requested modules to the list
  if cmdline.modules is not None:
    yaml_module_list=list(cmdline.modules)
    requestfile.register_classes(yaml_module_list)

  #Initialize a RequestFileRequest and run it
  req=requestfile.RequestFileRequest(requestfile=requestfile.fullpath)
  req.run()

