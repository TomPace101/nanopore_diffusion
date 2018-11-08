"""Run dolfin-convert"""

#This package
from ..requesthandler.shell import ShellCommandRequestBase
from ..requesthandler.yaml_manager import read as readyaml
from ..requesthandler import locators

#Locators
locators.folder_structure.update(mesh_xmlfile=['mesh',0,'xml'])

_DolfinConvertRunner_props_schema_yaml="""#DolfinConvertRunner
name: {type: string}
mshfile: {type: path}
xmlfile: {type: path}"""

class DolfinConvertRunner(ShellCommandRequestBase):
  """Run dolfin-convert
  
  User-defined attributes:
  
    - mshfile: Path to input .msh file
    - xmlfile = Path to output .xml file
       Two other files are also created by dolfin-convert, in the same directory,
       which contain mesh function data.
       This is the path to output file containing the mesh itself.
  """
  __slots__=('mshfile','xmlfile')
  _self_task=True
  _required_attrs=['mshfile','xmlfile']
  _inputfile_attrs=['mshfile']
  _outputfile_attrs=['xmlfile']
  _props_schema=readyaml(_DolfinConvertRunner_props_schema_yaml)
  @property
  def cmd_str(self):
    return "dolfin-convert %s %s"%(self.mshfile,self.xmlfile)