"""Pull together results from different requests or files into a single dataframe"""

#Standard library

#Site packages
import pandas as pd

#This package
from ..requesthandler.customization import CustomizableRequest, make_schema
from ..requesthandler import yaml_manager
from ..requesthandler import locators

#Locators
locators.folder_structure.update(CollectedFrame=['postproc',0])

_RawCollectionRequest_props_schema_yaml="""#RawCollectionRequest
outpath:
  type: pathlike
definitions:
  type: array
  items:
    type: array
    items:
      - {type: object}
      - type: array
        items: {type: pathlike}
calculations:
  type: array
  items:
    type: array
    items:
      - {type: string}
      - {type: object}
"""

class RawCollectionRequest(CustomizableRequest):
  """Collect data from an explicit list of files
  
  User-defined attributes:
  
    - outpath = path to store the output file containing the dataframe
    - definitions = sequence of pairs (mapping, file_list) where
      - mapping = dictionary {column name: dot-separated location of data within a loaded dictionary}
      - file_list = list of files to be loaded (one row per file, one loaded dictionary per file)
    - calculations = command sequence defining calculations to be done after the dataframe is constructed from the input files"""
    
  _self_task=True
  _required_attrs=['name','outpath','definitions']
  _outputfile_attrs=['outpath']
  _config_attrs=['definitions']
  _props_schema=make_schema(_RawCollectionRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(HDF5ConvertRequest, self).__init__(**kwargs)
    #Compile the input files
    self._more_inputfiles=[]
    for mapping,file_list in self.definitions:
      self._more_inputfiles += file_list
  def run(self):
    ##TODO!
    pass
  def do_calcs(self):
    ##TODO!
    pass

_CollectionRequest_props_schema_yaml="""#CollectionRequest
outpath:
  type: pathlike
requests:
  type: array
  items:
    type: request
definitions:
  type: array
  items:
    type: array
    items:
      - {type: object}
      - type: array
        items: {type: locator}
raw_definitions:
  type: array
  items:
    type: array
    items:
      - {type: object}
      - type: array
        items: {type: pathlike}
calculations:
  type: array
  items:
    type: array
    items:
      - {type: string}
      - {type: object}
"""

class CollectionRequest(CustomizableRequest):
  """Collect data from the output of other requests
  
  User-defined attributes:
  
    - outpath = path to store the output file containing the dataframe
    - requests = list of requests used to find the file to be read, as described below
    - definitions = sequence of pairs (mapping, locator_list) where
      - mapping = dictionary {column name: dot-separated location of data within a loaded dictionary}
      - locator_list = list of locators defining to be loaded (one row per file, one loaded dictionary per file)
    - calculations = command sequence defining calculations to be done after the dataframe is constructed from the input files
  
  This request finds files to include in the dataframe based on the list of requests, and the locator_list(s).
  For each locator in a given locator_list:
    For each request in the request list:
      Walk through all of its children, and their children, etc., and for each of those requests:
        If the request has ``_self_task`` set to True, the locator is rendered for that request.
  
  Calculated attributes:
  
    - raw_definitions = list of definition pairs suitable for RawCollectionRequest"""
    
  _self_task=True
  _required_attrs=['name','outpath','definitions','requests']
  _outputfile_attrs=['outpath']
  _config_attrs=['definitions','requests','calculations']
  _props_schema=make_schema(_CollectionRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(HDF5ConvertRequest, self).__init__(**kwargs)
    #Compile the raw definitions and input files
    self._more_inputfiles=[]
    self.raw_definitions=[]
    for mapping,locator_list in self.definitions:
      working_filelist=[]
      for req in self.requests:
        for ch in req.recursive_children():
          if getattr(ch,'_self_task',False):
            this_filelist = [ch.render(loc) for loc in locator_list]
            self._more_inputfiles += this_filelist
            working_filelist += this_filelist
      self.raw_definitions.append((mapping,working_filelist))
  def run(self):
      kwargs={'name':self.name+"_raw",'outpath':self.outpath,'definitions':self.raw_definitions}
      if getattr(self,'calculations',None) is not None:
        kwargs['calculations']=self.calculations
      RawCollectionRequest.complete(**kwargs)

#Register for loading from yaml
yaml_manager.register_classes([RawCollectionRequest, CollectionRequest])
