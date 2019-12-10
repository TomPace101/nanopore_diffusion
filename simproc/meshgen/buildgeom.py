"""Geneate gmsh .geo file(s) from standard geometric layout."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.request import make_schema, Request
from ..requesthandler import schema
from ..requesthandler import locators

#Locators
locators.folder_structure.update(geotemplate=['mesh','templates']) ##TODO

_GeometryDefinition_props_schema_yaml="""#GeometryDefinition
dimensions:
  type: integer
  minimum: 2
  maximum: 3
tmplfile: {type: pathlike}
tmplvars: {type: object}
outvars: {type: array}
ptdict: {type: object}
geomtable: {type: object}
surfloops: {type: object}
nonplanar: {type: array}
"""

class GeometryDefinition(schema.SelfValidating):
  """Geometry definition parameters
  
  Attributes:

    - dimensions = number of spatial dimensions (i.e. 2 for 2D mesh, 3 for 3D mesh)
    - tmplfile = geometry template file
    - tmplvars = mesh parameter variables needed by the geometry template file
    - outvars = list of gmsh variables to be output to the mesh metadata file (note: template variables are not always provided directly to gmsh)
    - ptdict = dictionary of points and their corresponding mesh density parameter name
    - geomtable = mapping of surfaces to sequence points
    - surfloops = mapping of surface loops to sequence of surfaces
    - nonplanar = list of surfaces that are not planar surfaces"""
  _required_attrs=['dimensions', 'tmplfile', 'tmplvars', 'outvars', 'ptdict', 'geomtable', 'surfloops', 'nonplanar']
  _props_schema=schema.SelfValidating.update_props_schema(_GeometryDefinition_props_schema_yaml)

schema.extra_types_dict['GeometryDefinition']=(GeometryDefinition,)


_BuildGeomRequest_props_schema_yaml="""#BuildGeomRequest
geomdef: {type: GeometryDefinition}
parameters: {type: object}
geofile: {type: pathlike}
"""

class BuildGeomRequest(Request):
  """Generate a gmsh .geo file from a specified geometric layout
  
  User-defined attributes:
  
    - geomdef = instance of GeometryDefinition
    - parameters = dictionary of parameter input values
    - geofile = path to output .geo file
  """
  _self_task=True
  _required_attrs=['name','geomdef','parameters','geofile']
  _outputfile_attrs=['geofile']
  _inputfile_attrs=[] ##TODO
  _config_attrs=['geomdef','parameters'] ##TODO: will this work as intended?
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
yaml_manager.register_classes([GeometryDefinition, BuildGeomRequest])
