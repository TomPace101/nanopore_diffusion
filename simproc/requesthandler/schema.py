"""Define class that can validate itself against a json schema"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility

#Site packages
import jsonschema

#This package
from . import filepath
from . import yaml_manager
from . import locators
from . import nested

#Validation partial setup (some setup must wait for Request class to be defined)
ValidatorClass = jsonschema.Draft4Validator
loc_class_tup=(locators.DataFile,locators.NameDelegator,locators.Delegator)
#jsonschema 2.6
extra_types_dict={'path':filepath.Path,
                  'locator':loc_class_tup,
                  'pathlike':(str,filepath.Path)+loc_class_tup, #This isn't the same thing as "pathlike" in python.org documentation
                  'array':(list,tuple),
                  'attrpath':(str,list,tuple),
                  'stored':nested.Stored}

def validation_error_string(err):
  "Return a string explaining the given validation error"
  #Get the basic message
  s=err.message
  #Provide a path if necessary
  if len(err.path)>0:
    s="%s: %s"%('.'.join([str(itm) for itm in err.path]),s)
  return s

#Basic schema for self-validation
#Note that validation is not applied to class attributes,
#but some of these could be instance attributes in some cases.
#Note also that properties inherited from nested.WithNested must also be included in the schema.
_SelfValidating_props_schema_yaml="""#SelfValidating
name:
  anyOf:
    - {type: 'null'}
    - {type: string}
store_globally: {type: boolean}
_required_attrs:
  type: array
  items: {type: string}
_props_schema: {type: object}
_additional_properties_ok: {type: boolean}
"""

class SelfValidating(nested.WithNested):
  """A class that can validate itself against a json schema
  
  Supported attributes:
  
    - _required_attrs: list of attribute names that must be defined for the request to be valid
    - _props_schema: jsonschema used to validate request configuration, as a dictionary
        The schema is for the 'properties' element only. The rest is provided internally.
    - _additional_properties_ok: boolean, True if the schema should allow additional properties"""
  _props_schema=yaml_manager.readstring(_SelfValidating_props_schema_yaml)
  _additional_properties_ok=False
  def __init__(self,**kwargs):
    #Validate kwargs
    if hasattr(self,'_props_schema'):
      self.validate_kwargs(**kwargs)
    #Load the attributes specified
    super(SelfValidating, self).__init__(**kwargs)
  @classmethod
  def update_props_schema(cls,yaml_str):
    """Return the property schema for a subclass
    
    The function is intended to be called by the subclass,
    possibly using super() to determine the appropriate base class
    (which, of course, this method will belong to).
    
    Arguments:
    
      - yaml_str = string containing yaml defining updates to _props_schema"""
    sub_schema=yaml_manager.readstring(yaml_str)
    schema={}
    schema.update(cls._props_schema)
    schema.update(sub_schema)
    return schema
  def _validation_schema(self):
    """Return the jsonschema validation schema for instances of this class"""
    return {'type':'object',
            'properties':self._props_schema,
            'required':getattr(self,'_required_attrs',[]),
            'additionalProperties':self._additional_properties_ok}
  def additional_validation(self,**kwargs):
    """Perform additional validation of the object data, beyond just the schema check

    The method provided in this base class is meant to be overridden by subclasses that need it.

    Arguments:

      - \*\*kwargs = arguments dictionary of attributes

    Returns:

      - list of error strings (including their indentation)"""
    return []
  def validate_kwargs(self,**kwargs):
    if hasattr(self,'_props_schema'):
      schema=self._validation_schema()
      validator=ValidatorClass(schema,types=extra_types_dict)
      errlist=["  - %s"%validation_error_string(err) for err in validator.iter_errors(kwargs)]
      errlist+=self.additional_validation(**kwargs)
      if len(errlist)>0:
        #Found errors: raise exception listing them all
        errlist.sort()
        errstr="Errors found in %s.\n"%type(self).__name__
        keylist=list(kwargs.keys())
        keylist.sort()
        errstr+='Received arguments:\n'
        errstr+='\n'.join(['  - %s: %s'%(k,kwargs[k]) for k in keylist])
        errstr+='\nErrors:\n'
        errstr+='\n'.join(errlist)
        raise Exception(errstr)
  def validate(self):
    d=self.to_dict(False)
    self.validate_kwargs(**d)
  