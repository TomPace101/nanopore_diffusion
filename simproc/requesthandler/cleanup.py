"For cleaning up the output files generated by other requests"

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os

#Site packages

#This package
from . import locators
from . import request
from . import yaml_manager

def cleanpath(target):
  """Remove the requested file or directory, and its parent if empty, recursively"""
  #Confirm that the path exists
  if target.exists():
    #File or folder
    if target.isFile:
      target.unlink()
    else:
      #Confirm that the directory is empty
      if len([p for p in target.iterdir()])==0:
        #Remove empty directory
        target.rmdir()
  #Check the parent directory, unless we are up high enough that we don't need to.
  #We do this even if the given path doesn't exist, because it might be that all the files were already deleted
  #and we only need to clean up the directories
  rp=target.relpath(locators.DATAFOLDER)
  min_len= 1 if rp.is_absolute() else 0
  if len(rp.parts) > min_len:
    cleanpath(target.parent)

_OutputCleanupRequest_props_schema_yaml="""#OutputCleanupRequest
clean: {type: array}
pathlist: {type: array}"""

class OutputCleanupRequest(request.Request):
  """Request to clean up the output of other requests
  
  User-Provided Attributes:
  
    - clean: a sequence of requests, whose output files are to be cleaned
    
  Calculated Attributes:
  
    - pathlist: list of files to be deleted"""
  _props_schema=request.make_schema(_OutputCleanupRequest_props_schema_yaml)
  _required_attrs=['clean']
  _self_task=True
  def __init__(self,**kwargs):
    #Initialization from base class
    super(OutputCleanupRequest, self).__init__(**kwargs)
    self.pathlist=[]
    #Iterate over provided list of requests
    for parent_req in self.clean:
      #Iterate over all the children of this request
      for child_req in parent_req.recursive_children():
        #If it creates a task, get the output files
        if child_req._self_task:
          self.pathlist.extend(child_req.outputfiles)
  def pre_run(self):
    """Final checks and preparatory steps"""
    #Confirm validation
    self.validate()
    #We don't need to assure output directories, as they might be deleted anyway
  def run(self):
    """Delete all the output files that exist, and remove empty directories"""
    #Final checks and preparatory steps
    self.pre_run()
    #Delete the listed paths
    for fpath in self.pathlist:
      cleanpath(fpath)

#Register for loading from yaml
yaml_manager.register_classes([OutputCleanupRequest])
