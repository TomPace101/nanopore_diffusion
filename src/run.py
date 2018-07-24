"""Command-line support for running requests"""

#Standard library
import argparse

#Site packages
import yaml
##from ruamel.yaml import YAML
##yaml=YAML()

#Local
import folderstructure as FS
import filepath

#Handle command-line execution
if __name__ == '__main__':
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('requestfile', help="Path to file containing the request(s) to run")
  #TODO: allow selecting a subset of the requests?
  cmdline=parser.parse_args()
  requestfile=filepath.Path(cmdline.requestfile,isFile=True)
  
  #Confirm that specified request file exists
  assert requestfile.exists(), "Could not find specified request file %s"%requestfile.fullpath
  
  #Load all of the requests
  ##TODO
  
