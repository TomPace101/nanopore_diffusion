#For generating parameter sets

#Standard library
import argparse
import itertools
import os
import os.path as osp
import sys

#Site packages
from jinja2 import Environment, FileSystemLoader

#Local
from folderstructure import *
import useful

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
      calcfields = sequence of field calculation specifiers, each specifier a pair:
        (calcfunc, kwargs)
        The calculation functions are also passed the fields dictionary as it exists at that point.
        Calculated fields are done after all other fields are known,
        and are executed in the order provided so that subsequent calculations can use them.
    To be generated by methods:
      tmpl_input = template input dictionary (see below)
  The template will be provided with the following input:
      allfields = a sequence of fields dictionaries,
        where each field"""
  __slots__=['outfile','tmplfile','constfields','rangefields','otherfiles','calcfields']

  def generate_fields(self):
    """Generator that calculates all the field dictionaries"""
    #Read in all the other files
    if getattr(self,'otherfiles',None) is not None:
      ##TODO
    #Loop through
    counter=1
    ##TODO: this is where it all really happens

  def prepare_tmpl_input(self):
    """Prepare the input needed by the template"""
    allfields=[field for field in self.generate_fields()]
    self.tmpl_input={'allfields':allfields}

  def do_generation(self):
    """Perform the requested parameter generation"""

    #Get the input dictionary for the template
    self.prepare_tmpl_input()

    #Load template
    env=Environment(loader=FileSystemLoader([pgtemplates_folder_folder,'.']), trim_blocks=True)
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
