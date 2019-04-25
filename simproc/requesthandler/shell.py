"""Support for requests that execute shell commands."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call

#This package
from . import request
from . import yaml_manager

class ShellCommandRequestBase(request.Request):
  """Base class for requests that execute shell commands
  
  Subclasses must do the following:
  
  - define _outputfile_attrs and/or _more_outputfiles, either as class or instance attributes
  - define a property attribute cmd_str that provides the shell command to be executed, as a string"""

  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Run the shell command
    call(self.cmd_str,shell=True)

_GeneralShellCommandRequest_props_schema_yaml="""#GeneralShellCommandRequest
name: {type: string}
outfile: {type: pathlike}
errfile:
  anyOf:
    - {type: pathlike}
    - {type: 'null'}
command: {type: string}"""

class GeneralShellCommandRequest(ShellCommandRequestBase):
  """Request for simple shell commands
  
  User-defined attributes:
  
    - command: string representing command to execute
    - outfile: Path to output file
    - errfile: Path to error output file, or None to redirect to `outfile`"""
  _self_task=True
  _required_attrs=['name','outfile','command']
  _config_attrs=['outfile','errfile','command']
  _props_schema=request.make_schema(_GeneralShellCommandRequest_props_schema_yaml)
  _outputfile_attrs=['outfile']
  @property
  def cmd_str(self):
    cmd="%s >'%s' "%(str(self.command),str(self.outfile))
    if getattr(self,'errfile',None) is None:
      cmd += "2>&1"
    else:
      cmd += "2>'%s'"%self.errfile
    return cmd
  
#Register for loading from yaml
yaml_manager.register_classes([GeneralShellCommandRequest])
