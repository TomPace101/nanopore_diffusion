"""Requests used only for debugging purposes"""

#Standard library
import time

#This package
from . import request
from . import yaml_manager
from . import shell

_DummyRequest_props_schema_yaml="""#DummyRequest
name: {type: string}
test:
  anyOf:
    - {type: string}
    - {type: number}
    - {type: pathlike}"""

class DummyRequest(request.Request):
  """A type of request used only for debugging and demonstration purposes
  
  User-defined attributes:
  
    - test: test data, which is printed when the request is run
  
  Here are some illustrations of basic request operations.
  
  >>> dr=DummyRequest(name='example',test='this_is_the_data')
  >>> dr.run()
  this_is_the_data
  >>> keylist=list(dr.task_definition.keys())
  >>> keylist.sort()
  >>> vlist=[dr.task_definition[k] for k in keylist]
  >>> list(zip(keylist,vlist)) # doctest: +ELLIPSIS +NORMALIZE_WHITESPACE
  [('actions', [(<bound method DummyRequest.run of <...DummyRequest object at 0x...>>,)]),
   ('file_dep', []),
   ('name', 'example'),
   ('targets', []),
   ('uptodate', [<doit.tools.config_changed object at 0x...>])]
  >>> invalid=DummyRequest(not_allowed=True)
  Traceback (most recent call last):
    ...
  Exception: Errors found in DummyRequest.
  Received arguments:
    - not_allowed: True
  Errors:
    - 'name' is a required property
    - 'test' is a required property
    - Additional properties are not allowed ('not_allowed' was unexpected)
  
  Be aware that Requests are not immutable,
  so it is possible to create an initially valid request,
  then modify it into an invalid one.
  
  >>> drq=DummyRequest(name='switcharoo',test='this_is_fine') #valid request
  >>> drq.test=None #None is not an allowed value for the 'test' attribute
  >>> drq.validate()
  Traceback (most recent call last):
    ...
  Exception: Errors found in DummyRequest.
  Received arguments:
    - name: switcharoo
    - test: None
  Errors:
    - test: None is not valid under any of the given schemas"""
  _self_task=True
  _config_attrs=['test']
  _required_attrs=['name','test']
  _props_schema=request.make_schema(_DummyRequest_props_schema_yaml)
  def run(self):
    self.validate()
    print(self.test)

_DummyShellRequest_props_schema_yaml="""#DummyShellRequest
name: {type: string}
outfile: {type: pathlike}
test:
  anyOf:
    - {type: string}
    - {type: number}
    - {type: pathlike}"""

class DummyShellRequest(shell.ShellCommandRequestBase):
  """A debugging test for a shell request
  
  User-defined attributes:
  
    - test: test data, which is passed to 'echo' when the request is run
    - outfile: Path to output file"""
  _self_task=True
  _required_attrs=['name','outfile','test']
  _config_attrs=_required_attrs
  _props_schema=request.make_schema(_DummyShellRequest_props_schema_yaml)
  _outputfile_attrs=['outfile']
  @property
  def cmd_str(self):
    return "echo '%s' >%s"%(str(self.test),self.renderstr(self.outfile))

class DummyShellAppendRequest(DummyShellRequest):
  @property
  def cmd_str(self):
    return "echo '%s' >>%s"%(str(self.test),self.renderstr(self.outfile))

_SleepRequest_props_schema_yaml="""#SleepRequest
name: {type: string}
delay:
  type: number
  minimum: 0
"""

class SleepRequest(request.Request):
  """Request to sleep for the specified number of seconds

  User-defined attributes:

    - delay: number of seconds to sleep, passed directly to time.sleep"""
  _self_task=True
  _required_attrs=['name','delay']
  _props_schema=request.make_schema(_SleepRequest_props_schema_yaml)
  def run(self):
    time.sleep(self.delay)

#Register for loading from yaml
yaml_manager.register_classes([DummyRequest, DummyShellRequest, DummyShellAppendRequest, SleepRequest])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
