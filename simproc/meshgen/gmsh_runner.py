"""Run gmsh"""

#This package
from ..requesthandler.shell import ShellCommandRequestBase
from ..requesthandler.yaml_manager import register_classes, read as readyaml
from ..requesthandler import locators

#Locators
locators.folder_structure.update(geotemplate=['mesh','templates'])
locators.folder_structure.update(geofile=['mesh',0,'geo'])
locators.folder_structure.update(mshfile=['mesh',0,'msh'])
locators.folder_structure.update(gmsh_outfile=['mesh',0,'gmsh_out'])
locators.folder_structure.update(meshmetafile=['mesh',0,'metadata'])

_GmshRequest_props_schema_yaml="""#GmshRequest
name: {type: string}
integer_arg: {type: integer}
geofile: {type: pathlike}
mshfile: {type: pathlike}
gmsh_outfile: {type: pathlike}
meshmetafile: {type: pathlike}"""

class GmshRequest(ShellCommandRequestBase):
  """Run gmsh
  
  User-defined attributes:
  
    - geofile: Path to input .geo file
    - mshfile: Path to output .msh file
    - gmsh_outfile: Path to text file to store gmsh message output
    - meshmetafile: optional, Path to yaml file to store mesh metadata
    - integer_arg: optional, integer to pass to gmsh on the command line, to specify meshing dimension
        defaults to 0, which indicates that the .geo file contains the appropriate `Mesh` command."""
  _self_task=True
  _required_attrs=['geofile','mshfile','gmsh_outfile']
  _outputfile_attrs=['mshfile','gmsh_outfile','meshmetafile']
  _inputfile_attrs=['geofile']
  _props_schema=readyaml(_GmshRequest_props_schema_yaml)
  @property
  def cmd_str(self):
    #Integer argument
    int_arg = getattr(self,'integer_arg',None)
    int_arg = 0 if int_arg is None else int_arg
    cmd = "gmsh -%d"%int_arg
    #meshmetafile, if provided
    if getattr(self,'meshmetafile',None) is not None:
      cmd += " -setstring meshmetafile '%s'"%self.meshmetafile
    #All the other files
    cmd += " -o '%s' '%s' >'%s'"%(self.mshfile,self.geofile,self.gmsh_outfile)
    return cmd

#Register for loading from yaml
register_classes([GmshRequest])
