"""Alter the base path of absolute paths in a doit database"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse

#Site packages


#Top-level function
def do_alteration(inputfile,outputfile,oldroot,newroot):
  pass

#Function to run from command line
def command_line_run():
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('inputfile',help="Path to input doit database file")
  parser.add_argument('outputfile',help="Path to output doit database file")
  parser.add_argument('oldroot',help="Root path within the input database file")
  parser.add_argument('newroot',help="Root path within the output database file")
  #Run
  do_alteration(parser.inputfile,parser.outputfile,parser.oldroot,parser.newroot)
  #Done
  return

#Handle command-line execution
if __name__ == '__main__':
  command_line_run()