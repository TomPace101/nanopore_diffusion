"""Support for requests that fill in templates."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call

#Site packages
from jinja2 import Environment

#This package
from . import yaml_manager
from . import customization

_TemplateFileRequest_props_schema_yaml="""#GeneralShellCommandRequest
name: {type: string}
tmplfile: {type: pathlike}
outfile: {type: pathlike}
data: {type: object}"""

class TemplateFileRequest(customization.CustomizableRequest):
  """General and base class for requests to fill in jinja2 template files
  
  You can subclass this by overriding the `get_template_input` method.
  You can also get the same effect using the customization interface
  to load a replacement `get_template_input` method from a module.
  Otherwise, the template input data is assumed to reside in `data`.
  
  User-defined attributes:
  
    - tmplfile: path to the input template file, as Path or string
    - outfile: path to the output file, as Path or string
    - data: dictionary of data used to compute the template input values"""
  _self_task=True
  _config_attrs=('tmplfile','outfile','data','modules','initializations','extra')
  _inputfile_attrs=['tmplfile']
  _outputfile_attrs=['outfile']
  _required_attrs=['name','tmplfile','outfile','data']
  _props_schema=customization.make_schema(_TemplateFileRequest_props_schema_yaml)
  def get_template_input(self):
    return self.data
  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Load the template
    with open(self.renderstr(self.tmplfile),'r') as fh:
      tdata=fh.read()
    env=Environment(extensions=['jinja2.ext.do'],trim_blocks=True,keep_trailing_newline=True)
    tmpl=env.from_string(tdata)
    #tmpl=Template(tdata,trim_blocks=True,keep_trailing_newline=True)
    #Do the calculations for the template values
    input_data=self.get_template_input()
    #Apply the data to the template
    out_data=tmpl.render(**input_data)
    #Write the output file
    with open(self.renderstr(self.outfile),'w') as fh:
      fh.write(out_data)

#Register for loading from yaml
yaml_manager.register_classes([TemplateFileRequest])
