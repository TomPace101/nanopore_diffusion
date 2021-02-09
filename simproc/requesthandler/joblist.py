"""Requests that generate a list of jobs to run"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import itertools

#Site packages
import pandas as pd

#This package
from . import yaml_manager
from . import locators
from . import commandseq
from . import logging

logger=logging.getLogger(__name__)

#Locators
locators.folder_structure.update(JobFile=['jobs'])

_JobListRequest_props_schema_yaml="""#JobListRequest
name: {type: string}
id_field: {type: string}
id_format: {type: string}
start_id: {type: integer}
outfile: {type: pathlike}
constfields: {type: object}
rangefields:
  type: object
  additionalProperties: {type: array}
prepcommands: {type: array}
postcommands: {type: array}
"""

class JobListRequest(commandseq.WithCommandsRequest):
  """Generate a list of jobs based on parametric variations

  User-defined attributes:

    - rangefields = optional dictionary of fields that will step through a sequence of values:
      {fieldname: [value, value, ...], ...}
    - otherlists = dictionary of other job lists to be used for creating variations:
      {joblist_attrpath: {old_fieldname: new_fieldname}, ...}


    - constfields = optional, dictionary of fields that will always be the same:
      {fieldname: value, ...}
    - outfile = optional, path to the output file, as Path or string
    - prepcommands = optional, sequence of commands to execute before template generation (e.g. to load additional data)
    - postcommands = optional, sequence of commands to execute after template generation (e.g. to output additional data)
    - id_field = optional, string for the job ID field name, defaults to "job_id"
    - id_format = optional, format for the job ID string, defaults to "%04d"
    - start_id = optional, starting job ID number, defaults to 1
  """
  _self_task=True
  _config_attrs=('constfields','rangefields','prepcommands','postcommands','outfile','id_format','start_id')
  _outputfile_attrs=['outfile']
  _validation_schema=commandseq.WithCommandsRequest.update_schema(_JobListRequest_props_schema_yaml)
  _validation_schema.required=[]
  def __init__(self,**kwargs):
    #Initialization from base class
    super(JobListRequest, self).__init__(**kwargs)
    #Get input and output files from the command sequences
    self.init_command_sequence('prepcommands')
    self.init_command_sequence('postcommands')
    #Apply default values
    self.id_field=getattr(self,"id_field","job_id")
    self.id_format=getattr(self,"id_format","%04d")
    self.start_id=getattr(self,"start_id",1)
    self.constfields=getattr(self,"constfields",{})
  def run(self):
    logger.debug("Running Request",request_class=type(self).__name__,request_name=getattr(self,"name",None))
    #Final checks and preparatory steps
    self.pre_run()
    #Run prepcommands
    self.process_command_sequence(attrpath='prepcommands',singlefunc=None,positional=False)
    #Process constant fields
    fieldnames_const=list(self.constfields.keys())
    values_const=list(self.constfields.values())
    #Set up iterator for the big loop (range fields)
    if getattr(self,'rangefields',None) is not None:
        fieldnames_range=tuple(self.rangefields.keys())
        iterator_range=itertools.product(*self.rangefields.values())
    else:
        #below, range_fields will need to be an empty dictionary
        fieldnames_range=tuple()
        iterator_range=[tuple()]
    #Construct the list of all field names (column headings) in order
    self.job_columns=[self.id_field]+fieldnames_const+list(fieldnames_range)
    #Initialize
    id_val=self.start_id
    self.joblist=[]
    #Loop over range fields    
    for values_range in iterator_range:
      #Initialize row
      job_id=self.id_format%id_val
      row=[job_id]+values_const+list(values_range)
      # range_fields=dict(zip(fieldnames_range,values_range))

      #Do calcfields
      # calclist=[tuple(*cf.items()) for cf in getattr(self,'calcfields',[])]
      result = True #needed for case of no calculations to be done
      # for funcname,kwargs in calclist:
      #   result = calcfuncs[funcname](fields,files_docs_dict,**kwargs)
      #   if not result:
      #     break
      #Check that document passes
      if result:
        self.joblist.append(row)
        id_val +=1


    #Dictionary for yaml output
    self.out_dict={'job_columns':self.job_columns,'jobs':self.joblist}
    #Convert job list to dataframe
    self.joblist_df=pd.DataFrame(self.joblist,columns=self.job_columns)
    #Output to yaml file if requested
    if getattr(self,'outfile',None) is not None:
      yaml_manager.writefile_flow(self.out_dict,self.render(self.outfile))
    #Run postcommands
    self.process_command_sequence(attrpath='postcommands',singlefunc=None,positional=False)
    #Done
    return

#Register for loading from yaml
yaml_manager.register_classes([JobListRequest])

