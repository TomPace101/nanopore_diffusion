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
dtype_outfile: {}
constfields: {type: object}
rangefields:
  type: object
  additionalProperties: {type: array}
otherlists:
  type: array
  items:
    type: array
    minItems: 2
    maxItems: 2
    items:
      - anyOf: #The attribute path
          - {type: string}
          - {type: array}
      - {type: object}
prepcommands: {type: array}
postcommands: {type: array}
"""

class JobListRequest(commandseq.WithCommandsRequest):
  """Generate a list of jobs based on parametric variations

  Note that none of the arguments below are required.
  But if none of them are provided, you'll just get any empty job list.

  User-defined attributes:

    - constfields = optional, dictionary of fields that will always be the same:
      {fieldname: value, ...}
    - rangefields = optional, dictionary of fields that will step through a sequence of values:
      {fieldname: [value, value, ...], ...}
    - otherlists = list of other job lists to be used for creating variations:
      [(joblist_attrpath, {old_fieldname: new_fieldname}), ...}
    - outfile = optional, path to the output file for the DataFrame, as Path or string
    - dtype_outfile = optional, path to the output file for the data types for the DataFrame, as Path or string
      (if you do use ``outfile``, this is strongly recommended as well)
    - prepcommands = optional, sequence of commands to execute before template generation (e.g. to load additional data)
    - postcommands = optional, sequence of commands to execute after template generation (e.g. to output additional data)
    - id_field = optional, string for the job ID field name, defaults to "job_id"
    - id_format = optional, format for the job ID string, defaults to "%04d"
    - start_id = optional, starting job ID number, defaults to 1

  The resulting job list is stored internally in the following calculated attributes:

    - job_columns: field names for the variables
    - joblist: list of data for each job
    - joblist_df: job list as a pandas DataFrame, constructed from self.joblist and self.job_columns
  """
  _self_task=True
  _config_attrs=('constfields','rangefields','otherlists','outfile','dtype_outfile'
                  'id_field','id_format','start_id','prepcommands','postcommands')
  _outputfile_attrs=['outfile','dtype_outfile']
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
    #Set up iterator for the range fields
    if getattr(self,'rangefields',None) is not None:
        fieldnames_range=list(self.rangefields.keys())
        def get_iterator_range():
          return itertools.product(*self.rangefields.values())
    else:
        #below, range_fields will need to be an empty dictionary
        fieldnames_range=[]
        def get_iterator_range():
          return [tuple()]
    #Process the list of other lists
    fieldnames_otherlists=[]
    fdict_seq=[]
    if getattr(self,"otherlists",None) is not None:
      iterators_otherlists_units=[]
      for jlpath,fieldsdict in self.otherlists:
        otherlist=self.get_nested(jlpath)
        nametupname=str(jlpath)
        fieldnames_otherlists+=list(fieldsdict.values())
        fdict_seq.append(fieldsdict)
        this_iterator=otherlist.itertuples(index=False,name=nametupname)
        iterators_otherlists_units.append(this_iterator)
      iterators_otherlists=itertools.product(*iterators_otherlists_units)
    else:
      iterators_otherlists=[tuple()]
    #Construct the list of all field names (column headings) in order
    self.job_columns=[self.id_field]+fieldnames_const+fieldnames_range+fieldnames_otherlists
    #Initialize
    id_val=self.start_id
    self.joblist=[]
    #Loop over fields from other job lists
    for source_tuples in iterators_otherlists:
      #Process the tuples from other lists
      otherlist_row_portion=[]
      if len(source_tuples)>0:
        for idx,ntup in enumerate(source_tuples):
          for old_name in fdict_seq[idx].keys():
            otherlist_row_portion.append(getattr(ntup,old_name))
      #Refresh the range iterator
      iterator_range = get_iterator_range()
      #Loop over range fields
      for values_range in iterator_range:
        #Initialize row
        job_id=self.id_format%id_val
        row=[job_id]+values_const+list(values_range)+otherlist_row_portion
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


    #Convert job list to dataframe
    self.joblist_df=pd.DataFrame(self.joblist,columns=self.job_columns)
    #Output to csv file if requested
    if getattr(self,'outfile',None) is not None:
      self.save_csv('joblist_df',self.outfile,index=False,dtype_csv_fpath=getattr(self,'dtype_outfile',None))
    #Run postcommands
    self.process_command_sequence(attrpath='postcommands',singlefunc=None,positional=False)
    #Done
    return

#Register for loading from yaml
yaml_manager.register_classes([JobListRequest])

