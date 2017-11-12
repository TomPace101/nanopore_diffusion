#For generating parameter sets

#Standard library
import argparse
from functools import reduce
import itertools
import operator
import os
import os.path as osp
import sys

#Site packages
from jinja2 import Environment, FileSystemLoader

#Local
from folderstructure import *
import useful

#Calculation functions
calcfuncs={}
def addcalcfunc(f):
  calcfuncs[f.__name__]=f

#Each calculation functiion must take fields and files_docs_dict as its positional arguments
#Reutrn False for fields that should not be added to the list of all fields
@addcalcfunc
def calc_lookup(fields,files_docs_dict,dest_field,src_field,ldict):
  """Look up value in a dictionary"""
  fields[dest_field]=ldict[fields[src_field]]
  return True

optable={'<':operator.lt, '<=':operator.le, '>':operator.gt, '>=':operator.ge,
         '==':operator.eq, '!=':operator.ne}
@addcalcfunc
def calc_comparison_ok(fields,files_docs_dict,opstr,field1,field2):
  """Do a comparison, and filter out this entry if it fails"""
  return optable[opstr](fields[field1],fields[field2])

@addcalcfunc
def calc_extremum(fields,files_docs_dict,dest_field,ismax,namelist):
  """Take max or min over fields specified by namelist"""
  vals=[fields[n] for n in namelist]
  if ismax:
    res=max(vals)
  else:
    res=min(vals)
  fields[dest_field]=res
  return True

@addcalcfunc
def calc_seq_split(fields,files_docs_dict,dest_field_seq,fieldname):
  """Split a field containing a sequence into multiple fields"""
  newfields=dict(zip(dest_field_seq,fields[fieldname]))
  fields.update(newfields)
  return True

#Get a value from a nested dictionary
def nested_location(obj,loc):
  """For nested dictionary obj, get item at location specified by loc
  obj={'a':{'b':{'c':11,'d':99}}}
  loc=['a','b','c']
  nested_location(obj,loc)
  11
  Note that loc must be a list, not a tuple."""
  return reduce(operator.getitem, [obj]+loc)

class ParameterGenerator(useful.ParameterSet):
  """ParameterSet used to generate other ParameterSets
  Attributes:
    To be read in:
      outfile = path (relative to params) for the parameter file to create
      tmplfile = path (relative to data/paramgen) for the template file
      constfields = dictionary of fields that will always be the same:
        {fieldname: value, ...}
      rangefields = dictionary of fields that will step through a sequence of values:
        {fieldname: [values, ...], ...}
      otherfiles = dictionary of other files, and mapping of field names from entries in that document:
        {filename: {fieldname: location in document, ...}, ...}
        Locations in the other document are specified by a sequence of names:
          [first name, second name, ...]
        File names must be paths relative to params.
        Empty field dictionaries may be used to provide documents for calcfields below.
      calcfields = sequence of field calculation specifiers, each specifier a dictionary:
        {calcfunc:kwargs}
        The calculation functions are also passed the fields dictionary as it exists at that point,
        and the complete yaml documents used to generate the otherfiles fields.
        Calculated fields are done after all other fields are known,
        and are executed in the order provided,
        so that the results of previous calculations can be used.
        The calculation functions can eliminate combinations by returning False.
    To be generated by methods:
      tmpl_input = template input dictionary (see below)
      num_generated = number of documents generated
  The template will be provided with the following input:
      allfields = a sequence of fields dictionaries"""
  __slots__=['outfile','tmplfile','constfields','rangefields','otherfiles','calcfields','tmpl_input','num_generated']

  def generate_fields(self):
    """Generator that calculates all the field dictionaries"""
    #Read in all the other files
    otherdocs={}
    if getattr(self,'otherfiles',None) is not None:
      for yfile in self.otherfiles.keys():
          otherdocs[yfile]=useful.readyaml_multidoc(osp.join(paramsfolder,yfile))
    #Set up iterator for the big loop
    if getattr(self,'rangefields',None) is not None:
        fieldnames_range=tuple(self.rangefields.keys())
        iterator_range=itertools.product(*self.rangefields.values())
    else:
        #range_fields needs to be an emtpy dictionary
        fieldnames_range=tuple()
        iterator_range=[tuple()]
    filenames_docs=tuple(otherdocs.keys())
    iterator_docs=itertools.product(*otherdocs.values())
    #Loop through
    counter=1
    for values_range in iterator_range:
      range_fields=dict(zip(fieldnames_range,values_range))
      for ydocs in iterator_docs:
        #Initialize the fields dictionary with a counter
        fields={'counter':counter}
        #Include the range fields
        fields.update(range_fields)
        #Include the constfields
        fields.update(getattr(self,'constfields',{}))
        #Get mapping of file to the document from that file
        files_docs_dict=dict(zip(filenames_docs,ydocs))
        if getattr(self,'otherfiles',None) is not None:
          #For each file
          for yfile,fieldmap in self.otherfiles.items():
            #Get the document to read from
            otherdoc=files_docs_dict[yfile]
            #For each destination field
            for fieldname,locseq in fieldmap.items():
              #Get the value at the specified location
              val=nested_location(otherdoc,locseq)
              #Put into specified field
              fields[fieldname]=val
        #Do calcfields
        calclist=[tuple(*cf.items()) for cf in getattr(self,'calcfields',[])]
        for funcname,kwargs in calclist:
          result = calcfuncs[funcname](fields,files_docs_dict,**kwargs)
          if not result:
            break
        #Check that document passes
        if result:
          yield fields
          counter +=1
    #Done at last
    self.num_generated = counter-1

  def prepare_tmpl_input(self):
    """Prepare the input needed by the template"""
    allfields=[fields for fields in self.generate_fields()]
    self.tmpl_input={'allfields':allfields}

  def do_generation(self):
    """Perform the requested parameter generation"""

    #Get the input dictionary for the template
    self.prepare_tmpl_input()

    #Load template
    env=Environment(loader=FileSystemLoader([pgtemplates_folder,'.']), trim_blocks=True)
    tmpl=env.get_template(self.tmplfile)

    #Render template
    doc = tmpl.render(self.tmpl_input)

    #Write output
    outfpath=osp.join(paramsfolder,self.outfile)
    os.makedirs(osp.dirname(outfpath),exist_ok=True)
    with open(outfpath,'w') as fp:
      fp.write(doc)

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  input_yaml_docstring="""File defining the parameter files to be generated
  This is a potentially multi-doc yaml file, where each document specifies one multi-doc yaml file to generate.
  Each document must provide the attributes for an instance of the ParameterGenerator class."""
  parser = argparse.ArgumentParser(description='Create parameter files')
  parser.add_argument('input_yaml', help=input_yaml_docstring)
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.input_yaml), "Input file does not exist: %s"%cmdline.input_yaml

  #Read in the yaml file
  alldocs=ParameterGenerator.all_from_yaml(cmdline.input_yaml)

  #Generate each requested parameter file
  for gendoc in alldocs:
    gendoc.do_generation()
    print("%s: %d documents"%(gendoc.outfile,gendoc.num_generated))
