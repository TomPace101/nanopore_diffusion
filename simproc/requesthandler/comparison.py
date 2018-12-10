"For comparing output to the expected output, to validate the code"

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import filecmp

#Site packages

#This package
from . import request
from . import yaml_manager

_FileComparisonRequest_props_schema_yaml="""#FileComparisonRequest
expected:
  anyOf:
    - {type: string}
    - {type: path}
received:
  anyOf:
    - {type: string}
    - {type: path}
"""

class FileComparisonRequest(request.Request):
  """Request to compare the contents of two files, and issue error if they don't match
  
  User-Provided Attributes:
  
    - expected: path to the file containing the expected output
    - received: path to the file containing the output produced by the code"""
  _props_schema=request.make_schema(_FileComparisonRequest_props_schema_yaml)
  _required_attrs=['expected','received']
  _inputfile_attrs=['expected','received']
  _self_task=True
  def run(self):
    args=(str(self.expected),str(self.received))
    ans=filecmp.cmp(*args,shallow=False)
    assert ans, "Found unexpected difference in files %s and %s"%args
    print("Files match: %s and %s"%args)

#Register for loading from yaml
yaml_manager.register_classes([FileComparisonRequest])
