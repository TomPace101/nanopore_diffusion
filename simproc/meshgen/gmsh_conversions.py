"""Run gmsh and subsequent conversion commands"""

#This package
from ..requesthandler.shell import ShellCommandRequestBase
from ..requesthandler.yaml_manager import read as readyaml

_GmshRunner_props_schema_yaml="""#GmshRunner
name: {type: string}
geofile: {type: path}
mshfile: {type: path}
txtfile: {type: path}
meshmetafile: {type: path}"""

class GmshRunner(ShellCommandRequestBase):
  """Run gmsh
  
  User-defined attributes:
  
    - geofile: Path to input .geo file
    - mshfile: Path to output .msh file
    - txtfile: Path to text file to store gmsh message output
    - meshmetafile: optional, Path to yaml file to store mesh metadata
    - integer_arg: optional, integer to pass to gmsh on the command line, to specify meshing dimension
        defaults to 0, which indicates that the .geo file contains the appropriate `Mesh` command."""
  __slots__=('geofile','mshfile','txtfile','meshmetafile','integer_arg')
  _self_task=True
  _required_attrs=['geofile','mshfile','txtfile']
  _outputfile_attrs=['mshfile','txtfile','meshmetafile']
  _inputfile_attrs=['geofile']
  _props_schema=readyaml(_GmshRunner_props_schema_yaml)
  _cmd_tmpl="gmsh -%d -setstring meshmetafile %s -o %s %s >%s"
  @property
  def cmd_str(self):
    int_arg=getattr(self,'integer_arg',None)
    int_arg = 0 if int_arg is None else int_arg
    ##TODO: meshmetafile is optional, so handle the case where it is missing
    ##TODO: do we need quotes around some file names?
    cmd=self._cmd_tmpl%(int_arg,self.meshmetafile,self.mshfile,self.geofile,self.txtfile)
    return cmd
