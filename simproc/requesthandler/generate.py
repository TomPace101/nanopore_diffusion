"""Support for requests that create their own children."""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import itertools

#Site packages

#This package
from . import yaml_manager
from .request import Request
from . import customization

_ParametricRequestListRequest_props_schema_yaml="""#ParametricRequestListRequest
request_type:
  anyOf:
    - {type: string}
    - {type: request}
constants: {type: object}
variations:
  type: object
  additionalProperties: {type: array}
other_parents:
  type: object
  additionalProperties: {type: request}
_children: {type: array}"""


class ParametricRequestListRequest(customization.CustomizableRequest):
  """General and base class for requests that parametrically generate child requests

  The process is:
  
    1) From the parent request attributes, construct the input fields for each child request
    2) From the input fields for each child request, compute the values of the keyword arguments to pass to its constructor
  
  Step 1 involves combining constant fields, fields that vary parametrically, and fields containing other requests.
  
  Step 2 involves calling the `get_child_kwargs` method.
  If the method returns None, then no request is generated for that combination of input fields.
  This can be used to filter out invalid combinations of parameters.

  You can subclass this by overriding the `get_child_kwargs` method.
  You can also get the same effect using the customization interface
  to load a replacement `get_child_kwargs` method from a module.
  Otherwise, the input fields are assumed to be identical to the child request keyword arguments.
  
  The fields passed to `get_child_kwargs` will also include 'index',
  which is the zero-based index number of the child in the sequence.
  
  User-defined attributes:
  
    - request_type = EITHER
          a type instance to use for the child requests
        OR
          a string containing the name of a request type registered with yaml_manager
    - constants = dictionary of variables that will be the same in each dictionary used to calculate request kwargs:
        {fieldname: value, ...}
    - variations = dictionary of fields that will vary over all possible combinations:
        {fieldname: [value, ...], ...}
        For example, {a: [1,2,3], b: [1,2,3]} will generate 9 combinations of a and b: (a=1,b=1), (a=1,b=2), (a=1,b=3), (a=2,b=1), ...
    - other_parents = dictionary of fields that will come from the children of other parent requests:
       {fieldname: parent_request, ...}
       For example, {other_request: parent} will create a field called `other_request`, which for each child request
       will contain one (immediate) child of the request called `parent`.
       Also, {request1: parent1, request2: parent2}, where parent1 and parent2 each have 3 child requests,
       will generate 9 combinations of request1 and request2.
       The iteration is only over the immediate children of the parent requests.
  
  Calculated Attributes:
  
    - _children = sequence of child requests"""
  _props_schema=customization.make_schema(_ParametricRequestListRequest_props_schema_yaml)
  _required_attrs=['request_type'] ##TODO
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
  def get_child_kwargs(self,index=None,**fields):
    """Compute a keyword arguments dictionary from the input dictionary"""
    outkwargs={}
    outkwargs.update(fields)
    return outkwargs
  def __init__(self,**kwargs):
    #Initialization from base class
    super(ParametricRequestListRequest, self).__init__(**kwargs)
    #Get a handle to the actual child request class
    if isinstance(self.request_type,str):
      childclass = yaml_manager.all_registered.get(self.request_type,None)
      assert childclass is not None, "Unable to find type %s in registered classes"%self.request_type
    else:
      assert isinstance(self.request_type,Request), "Invalid entry for request_type: %s of type %s."%(str(self.request_type),str(type(self.request_type)))
      childclass=self.request_type
    #Fields
    const_fields=getattr(self,'constants',{})
    variations=getattr(self,'variations',{})
    variation_fieldnames=list(variations.keys())
    variation_fieldnames.sort() #The order the fields were provided in might not have been preserved, so alphabetize.
    variation_fieldvalues=[variations[k] for k in variation_fieldnames]
    other_parents=getattr(self,'other_parents',{})
    op_fieldnames=tuple(other_parents.keys())
    op_c_generators=[p.all_children() for p in other_parents.values()] #The generator of all children for each other parent
    op_c_iterator=itertools.product(*op_c_generators) #Iterator over combinations of the other parents children
    #Loop through
    self._children=[]
    index=0
    for opc_tup in op_c_iterator:
      opc_fields=dict(zip(op_fieldnames,opc_tup))
      variation_iterator=itertools.product(*variation_fieldvalues)
      for variation_values in variation_iterator:
        variation_fields=dict(zip(variation_fieldnames,variation_values))
        #Put the fields together
        fields={'index':index}
        fields.update(const_fields)
        fields.update(variation_fields)
        fields.update(opc_fields)
        #Obtain arguments for child request constructor
        child_kwargs=self.get_child_kwargs(**fields)
        #Create child and add to list
        if child_kwargs is not None:
          ch_req=childclass(**child_kwargs)
          self._children.append(ch_req)
          index += 1
    return

#Register for loading from yaml
yaml_manager.register_classes([ParametricRequestListRequest])
