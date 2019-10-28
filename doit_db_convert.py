"""Alter the base path of absolute paths in a doit database"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import argparse
import dbm
import json

#Site packages

#Constants
ENCODING="utf-8"
DIRECT_COPY_KEYS=['checker:','_values_:']

def db_path_strip(infpath):
  return infpath[:-3]

def alter_task(indata,oldroot,newroot):
  "Do the alteration for a single task"
  outdata={}
  #Each key is processed differently
  for key,oldvalue in indata.items():
    if key in DIRECT_COPY_KEYS:
      #Direct copy
      outdata[key]=oldvalue
    elif key == 'deps:':
      #Each entry in the list is an absolute path
      newvalue=[fp.replace(oldroot,newroot) for fp in oldvalue]
      outdata[key]=newvalue
    else:
      #The key itself is an absolute path
      newkey=key.replace(oldroot,newroot)
      outdata[newkey]=oldvalue
  return outdata

#Top-level function
def do_alteration(inputfile,outputfile,oldroot,newroot):
  "Do the alteration for all the tasks in the database file"
  #Open the input and output files
  indb=dbm.open(db_path_strip(inputfile))
  outdb=dbm.open(db_path_strip(outputfile),'c')
  #Go through all the tasks
  for task_id in indb.keys():
    task_data=json.loads(indb[task_id].decode(ENCODING))
    new_task_data=alter_task(task_data,oldroot,newroot)
    outdb[task_id]=json.dumps(new_task_data).encode(ENCODING)
  #Close databases
  indb.close()
  outdb.close()
  #Done
  return

#Function to run from command line
def command_line_run():
  #Parse command line arguments
  parser = argparse.ArgumentParser(description=globals()['__doc__'])
  parser.add_argument('inputfile',help="Path to input doit database file")
  parser.add_argument('outputfile',help="Path to output doit database file")
  parser.add_argument('oldroot',help="Root path within the input database file")
  parser.add_argument('newroot',help="Root path within the output database file")
  cmdline=parser.parse_args()
  #Run
  do_alteration(cmdline.inputfile,cmdline.outputfile,cmdline.oldroot,cmdline.newroot)
  #Done
  return

#Handle command-line execution
if __name__ == '__main__':
  command_line_run()