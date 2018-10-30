"""Support for requests that fill in templates."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call

#Site packages
from jinja2 import Template

#This package
from . import request
from . import yaml_manager

_TemplateFileRequest_props_schema_yaml="""#GeneralShellCommandRequest
name: {type: string}
tmplfile: {type: path}
outfile: {type: path}
data: {type: object}"""


class TemplateFileRequest(request.Request):
  """General and base class for requests to fill in jinja2 template files
  
  You can subclass this by overriding the `get_template_input` method.
  Otherwise, the template input data is assumed to reside in `data`.
  
  User-defined attributes:
  
    - tmplfile: Path to the input template file
    - outfile: Path to the output file
    - data: dictionary of data used to compute the template input values"""
  __slots__=('tmplfile','outfile','data')
  _self_task=True
  _inputfile_attrs=['tmplfile']
  _outputfile_attrs=['outfile']
  _required_attrs=['tmplfile','outfile','data']
  _props_schema=yaml_manager.read(_TemplateFileRequest_props_schema_yaml)
  def get_template_input(self):
    return self.data
  def run(self):
    #Confirm validation
    self.validate()
    #Create directories for output files if necessary
    self.assure_output_dirs()
    #Load the template
    with open(self.tmplfile,'r') as fh:
      tdata=fh.read()
    tmpl=Template(tdata,trim_blocks=True)
    #Do the calculations for the template values
    input_data=self.get_template_input()
    #Write the output file
    ##TODO
