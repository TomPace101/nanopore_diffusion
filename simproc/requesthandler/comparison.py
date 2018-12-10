"""For comparing output to the expected output, to validate the code

Setup data for test

>>> files=['/tmp/a.txt', '/tmp/b.txt', '/tmp/c.txt']
>>> files_data=['hello world\\n']*2 + ['goodbye!\\n']
>>> for fpath,data in zip(files,files_data):
...   with open(fpath,'w') as fh:
...     fh.write(data)
...
12
12
9

Example with files that do match

>>> req1=FileComparisonRequest(expected=files[0],received=files[1])
>>> req1.run()
Files match: /tmp/a.txt and /tmp/b.txt

Example with files that don't match

>>> req2=FileComparisonRequest(expected=files[1],received=files[2])
>>> req2.run()
Traceback (most recent call last):
  ...
AssertionError: Found unexpected difference in files /tmp/b.txt and /tmp/c.txt

Clean up

>>> import os
>>> for fpath in files:
...   os.remove(fpath)
...

"""

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
