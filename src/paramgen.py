#For generating parameter sets

#Standard library
import argparse
import os
import os.path as osp
import sys

#Site packages

#Local
from folderstructure import *
import useful

class ParameterGenerator(useful.ParameterSet):
  """ParameterSet used to generate other ParameterSets
  Attributes:
    outfile = path (relative to params) for the parameter file to create"""
  __slots__=['outfile']

  def do_generation(self):
  """Perform the requested parameter generation"""

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
