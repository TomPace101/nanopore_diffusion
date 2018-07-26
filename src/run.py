"""Command-line support for running requests"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import importlib
import sys

#Site packages
from ruamel.yaml import YAML
yaml=YAML(typ="safe", pure=True)

#Local
import folderstructure as FS
import filepath
import request

#Complete list of all modules defining classes we want to load from yaml
yaml_module_list=['locators','request']

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
  
  #Initialize a RequestFileRequest and run it
  req=request.RequestFileRequest(requestfile=requestfile.fullpath)
  req.run()

  ##TODO: the stuff below probably goes into request.RequestFileRequest
  ##Except that request can't load other modules defining requests.
  ##So, instead, those need to move to a new module.

  
  #Add the requested modules to the list
  if cmdline.modules is not None:
    yaml_module_list+=list(cmdline.modules)
  
  #Load and register all classes that may be stored in yaml
  for modname in yaml_module_list:
    loaded_module=importlib.import_module(modname)
    for yclass in getattr(loaded_module,'yaml_classes',[]):
      yaml.register_class(yclass)
  
  #Get an iterator for the individual requests
  with open(requestfile.fullpath,'r') as fp:
    dat=fp.read()
  allreqs=yaml.load_all(dat)
  
  #Process each request
  for req in allreqs:
    req.run()
  
