"""Base class for requests that execute shell commands."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call

#This package
from . import request
from . import yaml_manager

class ShellCommandRequest(request.Request):
  """Base class for requests that execute shell commands
  
  Subclasses must do the following:
  
  - define _outputfile_attrs and/or _more_outputfiles, either as class or instance attributes
  - define a property attribute cmd_str that provides the shell command to be executed, as a string"""

  def run(self):
    #Confirm validation
    self.validate()
    #Create directories for output files if necessary
    allpaths=[getattr(self,oattr) for oattr in getattr(self,_outputfile_attrs,[])]
    allpaths+=getattr(self,_more_outputfiles,[])
    for fpath in allpaths:
      fpath.assure_dir()
    #Run the shell command
    call(self.cmd_str,shell=True)

