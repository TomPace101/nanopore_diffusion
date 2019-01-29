"""Run dolfin-convert"""

#This package
from ..requesthandler.shell import ShellCommandRequestBase
from ..requesthandler.yaml_manager import register_classes, read as readyaml
from ..requesthandler import locators

#Locators
locators.folder_structure.update(mesh_xmlfile=['mesh',0,'xml'])
locators.folder_structure.update(dconv_outfile=['mesh',0,'dconv_out'])

_DolfinConvertRequest_props_schema_yaml="""#DolfinConvertRequest
name: {type: string}
mshfile: {type: pathlike}
xmlfile: {type: pathlike}
txtfile: {type: pathlike}"""

class DolfinConvertRequest(ShellCommandRequestBase):
  """Run dolfin-convert
  
  User-defined attributes:
  
    - mshfile: Path to input .msh file
    - xmlfile = Path to output .xml file
       Two other files are also created by dolfin-convert, in the same directory,
       which contain mesh function data.
       This is the path to output file containing the mesh itself.
    - txtfile: Path to text file to store dolfin-convert message output
  """
  _self_task=True
  _required_attrs=['mshfile','xmlfile','txtfile']
  _inputfile_attrs=['mshfile']
  _outputfile_attrs=['xmlfile','txtfile']
  _props_schema=readyaml(_DolfinConvertRequest_props_schema_yaml)
  @property
  def cmd_str(self):
    return "dolfin-convert %s %s > %s"%(self.mshfile,self.xmlfile,self.txtfile)

#Register for loading from yaml
register_classes([DolfinConvertRequest])
