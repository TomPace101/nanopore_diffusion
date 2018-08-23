"""Run gmsh and subsequent conversion commands"""

#This package
from ..requesthandler.shell import ShellCommandRequest
from ..requesthandler.yaml_manager import read as readyaml

_GmshRunner_props_schema_yaml="""#GmshRunner
name: {type: string}
geofile: {type: path}
mshfile: {type: path}
txtfile: {type: path}
meshmetafile: {type: path}"""

class GmshRunner(ShellCommandRequest):
  """Run gmsh
  
  User-defined attributes:
  
    - geofile: Path to input .geo file
    - mshfile: Path to output .msh file
    - txtfile = Path to text file to store gmsh message output
    - meshmetafile = optional, Path to yaml file to store mesh metadata"""
  __slots__=('geofile','mshfile','txtfile','meshmetafile')
  _self_task=True
  _required_attrs=['geofile','mshfile','txtfile']
  _props_schema=readyaml(_GmshRunner_props_schema_yaml)
  
