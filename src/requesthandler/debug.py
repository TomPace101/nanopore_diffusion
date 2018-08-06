"""Requests used only for debugging purposes"""

#This package
from . import request
from . import yaml_manager

_DummyRequest_props_schema_yaml="""#DummyRequest
name: {type: string}
test:
  anyOf:
    - {type: string}
    - {type: number}
    - {type: path}"""

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
  [('actions', [(<bound method DummyRequest.run of <request.DummyRequest object at 0x...>>,)]),
   ('file_dep', []),
   ('name', 'example'),
   ('targets', []),
   ('uptodate', [<doit.tools.config_changed object at 0x...>])]
  >>> invalid=DummyRequest(not_allowed=True)
  Traceback (most recent call last):
    ...
  Exception: Errors found in DummyRequest:
    - 'name' is a required property
  <BLANKLINE>
  Failed validating 'required' in schema:
      {'additionalProperties': False,
       'properties': {'name': {'type': 'string'},
                      'test': {'anyOf': [{'type': 'string'},
                                         {'type': 'number'}]}},
       'required': ['name', 'test'],
       'type': 'object'}
  <BLANKLINE>
  On instance:
      {'not_allowed': True}
    - 'test' is a required property
  <BLANKLINE>
  Failed validating 'required' in schema:
      {'additionalProperties': False,
       'properties': {'name': {'type': 'string'},
                      'test': {'anyOf': [{'type': 'string'},
                                         {'type': 'number'}]}},
       'required': ['name', 'test'],
       'type': 'object'}
  <BLANKLINE>
  On instance:
      {'not_allowed': True}
    - Additional properties are not allowed ('not_allowed' was unexpected)
  <BLANKLINE>
  Failed validating 'additionalProperties' in schema:
      {'additionalProperties': False,
       'properties': {'name': {'type': 'string'},
                      'test': {'anyOf': [{'type': 'string'},
                                         {'type': 'number'}]}},
       'required': ['name', 'test'],
       'type': 'object'}
  <BLANKLINE>
  On instance:
      {'not_allowed': True}"""
  __slots__=('test')
  _self_task=True
  _config_attrs=['test']
  _required_attrs=['name','test']
  _props_schema=yaml_manager.read(_DummyRequest_props_schema_yaml)
  def run(self):
    print(self.test)

#Register for loading from yaml
yaml_manager.register_classes([DummyRequest])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
