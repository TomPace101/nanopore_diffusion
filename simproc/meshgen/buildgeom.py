"""Geneate gmsh .geo file(s) from standard geometric layout."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.request import make_schema, Request
from ..requesthandler import locators

#Locators
locators.folder_structure.update(geotemplate=['mesh','templates'])

_BuildGeomRequest_props_schema_yaml="""#BuildGeomRequest
"""

class BuildGeomRequest(Request):
  """Generate a gmsh .geo file from a specified geometric layout
  
  User-defined attributes:
  
    - 
  """
  _self_task=True
  _required_attrs=['name'] ##TODO
  _outputfile_attrs=[] ##TODO
  _inputfile_attrs=[] ##TODO
  _config_attrs=[] ##TODO
  _props_schema=make_schema(_BuildGeomRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(BuildGeomRequest, self).__init__(**kwargs)
    ##TODO
  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    ##TODO
    #Done
    return

#Register for loading from yaml
yaml_manager.register_classes([BuildGeomRequest])
