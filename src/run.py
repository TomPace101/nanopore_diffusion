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
  parser.add_argument('requestfile',nargs="+",help="Path to file containing the request(s) to run. Multiple request files may be listed.")
  parser.add_argument('--modules',nargs="+",metavar="MODULE",help="Additional python modules defining classes loadable from yaml input")
  #TODO: allow selecting a subset of the requests?
  cmdline=parser.parse_args()
  requestfiles=[filepath.Path(rf,isFile=True) for rf in cmdline.requestfile]
  
  #Confirm that specified request file(s) exist(s)
  nonexist=[rf for rf in requestfiles if not rf.exists()]
  assert len(nonexist)==0, "Could not find specified request file(s) %s"%str([rf.fullpath for fr in nonexist])

  #Add the requested modules to the list
  if cmdline.modules is not None:
    yaml_module_list=list(cmdline.modules)
    requestfile.register_classes(yaml_module_list)

  #Initialize a RequestFileListRequest and run it
  req=requestfile.RequestFileListRequest(requestfiles=[rf.fullpath for rf in requestfiles])
  req.run()

