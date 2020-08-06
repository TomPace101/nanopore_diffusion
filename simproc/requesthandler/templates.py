"""Support for requests that fill in templates."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from subprocess import call

#Site packages
from jinja2 import Environment, FileSystemLoader

#This package
from . import yaml_manager
from . import commandseq
from . import logging

logger=logging.getLogger(__name__)

_TemplateFileRequest_props_schema_yaml="""#TemplateFileRequest
name: {type: string}
searchpaths:
  type: array
  items: {type: pathlike}
tmplfile: {type: pathlike}
outfile: {type: pathlike}
data: {type: object}
prepcommands: {type: array}
postcommands: {type: array}
"""

class TemplateFileRequest(commandseq.WithCommandsRequest):
  """General and base class for requests to fill in jinja2 template files
  
  You can subclass this by overriding the `get_template_input` method.
  You can also get the same effect using the customization interface
  to load a replacement `get_template_input` method from a module.
  Otherwise, the template input data is assumed to reside in `data`.
  
  User-defined attributes:
  
    - searchpaths: optional list of folders to add to the jinja2 search path for template inheritance
    - tmplfile: path to the input template file, as Path or string
    - outfile: path to the output file, as Path or string
    - data: dictionary of data used to compute the template input values
    - prepcommands = sequence of commands to execute before template generation (e.g. to load additional data)
    - postcommands = sequence of commands to execute after template generation (e.g. to output additional data)"""
  _self_task=True
  _config_attrs=('tmplfile','outfile','data','modules','initializations','extra')
  _inputfile_attrs=['tmplfile']
  _outputfile_attrs=['outfile']
  _validation_schema=commandseq.WithCommandsRequest.update_schema(_TemplateFileRequest_props_schema_yaml)
  _validation_schema.required=['name','tmplfile','outfile','data']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(TemplateFileRequest, self).__init__(**kwargs)
    #Default command attributes
    if not hasattr(self,'prepcommands'):
      self.prepcommands=[]
    if not hasattr(self,'postcommands'):
      self.postcommands=[]
    #Get input files
    self._more_inputfiles=getattr(self,'_more_inputfiles',[]) #Initialize attribute if it doesn't already exist
    self._more_inputfiles+=[fp for k,fp in getattr(self,'loadfiles',{}).items()]
    self._more_inputfiles+=self.list_iofiles(self.prepcommands,['filename','infpath'],'_inputfiles')
    self._more_inputfiles+=self.list_iofiles(self.postcommands,['filename','infpath'],'_inputfiles')
    #Get output files
    self._more_outputfiles=getattr(self,'_more_outputfiles',[]) #Initialize attribute if it doesn't already exist
    self._more_outputfiles+=[figprops.outfpath for figprops in getattr(self,'figures',[])]
    self._more_outputfiles+=self.list_iofiles(self.prepcommands,['filename','outfpath'],'_outputfiles')
    self._more_outputfiles+=self.list_iofiles(self.postcommands,['filename','outfpath'],'_outputfiles')
  def get_template_input(self):
    return self.data
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Default search paths
    searchpaths=getattr(self,'searchpaths',[])
    #Load the template
    with open(self.renderstr(self.tmplfile),'r') as fh:
      tdata=fh.read()
    env=Environment(loader=FileSystemLoader(searchpaths),extensions=['jinja2.ext.do'],trim_blocks=True,keep_trailing_newline=True)
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
