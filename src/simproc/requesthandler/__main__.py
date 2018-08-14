"""Command-line support for running requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import doctest
import importlib
import sys

#Site packages
import pkg_resources #part of setuptools

#Local
from . import *

#Paths to files containing doctests
tutorial_file=pkg_resources.resource_filename(__name__,'tutorial.rst')

#Handle command-line execution
if __name__ == '__main__':
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('requestfile',nargs="*",help="Path to file containing the request(s) to run. Multiple request files may be listed.")
  parser.add_argument('--verbose',action='store_true',help="Provide verbose output where appropriate.")
  parser.add_argument('--modules',nargs="+",metavar="MODULE",help="Additional python modules defining classes loadable from yaml input")
  parser.add_argument('--validate',action='store_true',help="Perform validation. If requestfiles are also listed, validation is run first.")
  #TODO: allow selecting a subset of the requests?
  cmdline=parser.parse_args()
  
  #run validation if requested
  if cmdline.validate:
    reslist=[]
    reslist.append(doctest.testmod(filepath,verbose=cmdline.verbose))
    print("---")
    reslist.append(doctest.testmod(debug,verbose=cmdline.verbose))
    print("---")
    reslist.append(doctest.testfile(tutorial_file,module_relative=False,verbose=cmdline.verbose)) #To make sure the file can be found when the package is zipped, we have already found its absolute path
    print("---")
    fails,atts=[sum(l) for l in zip(*reslist)]
    print("Passed %d/%d total"%(atts-fails,atts))

  #Confirm that specified request file(s) exist(s)
  file_list=[filepath.Path(rf,isFile=True) for rf in cmdline.requestfile]
  nonexist=[rf for rf in file_list if not rf.exists()]
  assert len(nonexist)==0, "Could not find specified request file(s) %s"%str([rf.fullpath for rf in nonexist])

  #Load the requested modules
  if cmdline.modules is not None:
    customization.load_modules(cmdline.modules)

  #Initialize a RequestFileListRequest and run it
  req=requestfile.RequestFileListRequest(requestfiles=file_list)
  req.run()

