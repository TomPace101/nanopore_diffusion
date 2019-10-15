"""Pull together results from different requests or files into a single dataframe"""

#Standard library

#Site packages
import pandas as pd

#This package
from ..requesthandler.commandseq import WithCommandsRequest, make_schema
from ..requesthandler import yaml_manager
from ..requesthandler import pickle_manager
from ..requesthandler import locators
from ..requesthandler import nested

#Locators
locators.folder_structure.update(CollectedFrame=['postproc',0])

_RawCollectionRequest_props_schema_yaml="""#RawCollectionRequest
outpath:
  type: pathlike
multidoc: {type: boolean}
definitions:
  type: array
  items:
    type: object
    properties:
      mapping: {type: object}
      file_list:
        type: array
        items: {type: pathlike}
calculations:
  type: array
  items:
    type: array
    items:
      - {type: string}
      - {type: object}
"""

class RawCollectionRequest(WithCommandsRequest):
  """Collect data from an explicit list of files
  
  User-defined attributes:
  
    - outpath = path to store the output file containing the dataframe
    - multidoc = boolean, True if yaml files contain multiple documents, False if just one document each
    - definitions = sequence of dictionaries:
      - mapping = dictionary {column name: dot-separated location of data within a loaded dictionary} (each loaded dictionary comes from one file)
      - file_list = list of files to be loaded (one row per file)
      Note that every file_list must have the same number of files.
      Very often, you will have only a single definition, meaning each row gets all of its data from a single file.
      If each row needs to combine data from multiple files, you'll need more than one definition in the list.
      If multiple definitions include the same column, the order of the definitions will matter: later ones will overwrite earlier ones.
    - calculations = command sequence defining calculations to be done after the dataframe is constructed from the input files"""
    
  _self_task=True
  _required_attrs=['name','outpath','definitions']
  _outputfile_attrs=['outpath']
  _config_attrs=['definitions','multidoc']
  _props_schema=make_schema(_RawCollectionRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(RawCollectionRequest, self).__init__(**kwargs)
    #Compile the input files
    self._more_inputfiles=[]
    for mapping,file_list in self.definitions:
      self._more_inputfiles += [self.renderstr(fp) for fp in file_list]
    #Default multidoc
    if not hasattr(self,'multidoc'):
      self.multidoc=False
  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Get the number of rows
    nrow_list=[]
    for defn in self.definitions:
      nrow_list.append(len(defn['file_list']))
    nrows=nrow_list[0]
    assert all([nr==nrows for nr in nrow_list[1:]]), "Differing numbers of rows: %s"%str(nrow_list)
    #Initialize dictionaries to store data
    rows=[dict() for i in range(nrows)]
    #Populate dictionaries
    for defn in self.definitions:
      mapping=defn['mapping']
      file_list=defn['file_list']
      for fn,fpath in enumerate(file_list):
        #Load the file
        obj=yaml_manager.readfile(self.renderstr(fpath),self.multidoc)
        #Get the dictionary for the row
        flatdict=dict([(col,nested.get_nested(obj,dpath)) for col,dpath in mapping.items()])
        #Store
        rows[fn].update(flatdict)
    #Create dataframe
    self.df=pd.DataFrame(rows)
    #Do post-collection calculations
    self.process_command_sequence(attrpath='calculations',singlefunc=None,positional=False)
    #Save results
    pickle_manager.writefile(self.df,self.renderstr(self.outpath))

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
    type: object
    properties:
      mapping: {type: object}
      locator_list:
        type: array
        items: {type: locator}
raw_definitions:
  type: array
  items:
    type: object
    properties:
      mapping: {type: object}
      file_list:
        type: array
        items: {type: pathlike}
calculations:
  type: array
  items:
    type: array
    items:
      - {type: string}
      - {type: object}
"""

class CollectionRequest(WithCommandsRequest):
  """Collect data from the output of other requests
  
  User-defined attributes:
  
    - outpath = path to store the output file containing the dataframe
    - requests = list of requests used to find the file to be read, as described below
    - definitions = sequence of dictionaries:
      - mapping = dictionary {column name: dot-separated location of data within a loaded dictionary} (each loaded dictionary comes from one locator)
      - locator_list = list of locators for files to be loaded (one row per locator)
      Note that every locator_list must have the same number of locators.
      Very often, you will have only a single definition, meaning each row gets all of its data from a single file.
      If each row needs to combine data from multiple files, you'll need more than one definition in the list.
      If multiple definitions include the same column, the order of the definitions will matter: later ones will overwrite earlier ones.
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
  _config_attrs=['definitions','_more_inputfiles','calculations']
  _props_schema=make_schema(_CollectionRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    ##TODO: should this perhaps just use RawCollectionRequest as its children instead of calling `complete` on them when run?
    #Initialization from base class
    super(CollectionRequest, self).__init__(**kwargs)
    #Compile the raw definitions and input files
    self._more_inputfiles=[]
    self.raw_definitions=[]
    for defn in self.definitions:
      working_filelist=[]
      for req in self.requests:
        for ch in req.recursive_children():
          if getattr(ch,'_self_task',False):
            this_filelist = [ch.render(loc) for loc in defn['locator_list']]
            self._more_inputfiles += this_filelist
            working_filelist += this_filelist
      self.raw_definitions.append({'mapping':defn['mapping'],'file_list':working_filelist})
  def run(self):
      kwargs={'name':self.name+"_raw",'outpath':self.render(self.outpath),'definitions':self.raw_definitions}
      for argname in ['calculations','modules','methods']:
        if getattr(self,argname,None) is not None:
          kwargs[argname]=getattr(self,argname)
      RawCollectionRequest.complete(**kwargs)

_SimpleCollectionRequest_props_schema_yaml="""#SimpleCollectionRequest
outpath:
  type: pathlike
requests:
  type: array
  items:
    type: request
mapping: {type: object}
locator_list:
  type: array
  items: {type: locator}
calculations:
  type: array
  items:
    type: array
    items:
      - {type: string}
      - {type: object}
child: {type: object} ##actually a request, but to_dict is recursive, so the validation will receive a dictionary
"""

class SimpleCollectionRequest(WithCommandsRequest):
  """Collect data from the output of other requests
  
  User-defined attributes:
  
    - outpath = path to store the output file containing the dataframe
    - requests = list of requests used to find the file to be read, as described below
    - mapping = dictionary {column name: dot-separated location of data within a loaded dictionary} (each loaded dictionary comes from one locator)
    - locator_list = list of locators for files to be loaded (one table row per locator)
    - calculations = command sequence defining calculations to be done after the dataframe is constructed from the input files
  
  The generated table will use exactly one file to get the data in each table row.
  
  Calculated attributes:
  
    - child = the CollectionRequest created from this request"""
  _self_task=False
  _required_attrs=['name','outpath','requests','mapping','locator_list']
  _child_attrs=['child']
  _props_schema=make_schema(_SimpleCollectionRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(SimpleCollectionRequest, self).__init__(**kwargs)
    #Set up the child
    kwargs={
      'name':self.name+".child",
      'outpath':self.outpath,
      'requests':self.requests,
      'definitions':[{'mapping':self.mapping,'locator_list':self.locator_list}]
      }
    for argname in ['calculations','modules','methods']:
      if hasattr(self,argname):
        kwargs[argname]=getattr(self,argname)
    self.child=CollectionRequest(**kwargs)

#Register for loading from yaml
yaml_manager.register_classes([RawCollectionRequest, CollectionRequest, SimpleCollectionRequest])
