"For running requests simultaneously"

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import math
import os
import subprocess
import sys
import time

#This package
from . import request
from . import yaml_manager
from .filepath import Path

main_mod_name=request.__package__.split('.')[0]
modpath_parts=Path(__file__).parts
idx=modpath_parts.index(main_mod_name)
work_path=str(Path(*modpath_parts[:idx])) #Path to the parent directory of the top-level package

_SimultaneousRequestQueue_props_schema_yaml="""#SimultaneousRequestQueue
name: {type: string}
num_workers:
  type: integer
  minimum: 1
queue: {type: array}
delay:
  type: number
  minimum: 0
tmploc: {type: pathlike}
tmpfmt: {type: string}
"""

class SimultaneousRequestQueue(request.Request):
  """Run a queue of requests, with more than one allowed to run simultaneously
  
  User-defined attributes:
  
    - num_workers: number of requests to run in parallel
    - queue: sequence of the requests to run
    - delay: optional, time (in seconds), to wait between checks for request completion
    - tmploc: optional, path to folder to write temporary input files to
      Note that, at present, cleanup won't remove this folder if it is created. Sorry.
    - tmpfmt: optional, ugly, and I don't recommend you use it without looking at the source code"""
  _self_task=False
  _required_attrs=['num_workers','queue']
  _props_schema=request.make_schema(_SimultaneousRequestQueue_props_schema_yaml)
  _child_seq_attrs=['queue']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(SimultaneousRequestQueue, self).__init__(**kwargs)
    #Default values
    self.delay = getattr(self,'delay',30)
    self.tmploc = Path(getattr(self,'tmploc','.'),isFile=False).expanduser().resolve()
    self.tmpfmt = getattr(self,'tmpfmt',"req%0{}d.yaml")
  def run(self):
    #Make sure the temporary files directory exists
    self.tmploc.assure_dir()
    #Get the template for the temporary file names
    tmp_tmpl=self.tmpfmt.format(1+math.floor(math.log10(len(self.queue))))
    #Initialize
    running=[] #To store pairs (Popen, fpath)
    indx=0
    #Loop until the queue is completed
    while indx<len(self.queue) or len(running)>0:
      #Start new processes to keep the requested number of workers going simultaneously
      while len(running)<self.num_workers and indx<len(self.queue):
        #Take the next request off the queue
        req=self.queue[indx]
        #Write the input file
        fpath=self.tmploc / (tmp_tmpl%indx)
        yaml_manager.writefile([req],fpath)
        #Start the subprocess
        args=(sys.executable,'-m',main_mod_name,str(fpath))
        p=subprocess.Popen(args,cwd=work_path,shell=False)
        #Add to list of running processes
        running.append((p,fpath))
        #Prepare for next item from queue
        indx+=1
      #Wait until we check again
      if len(running)>0:
        time.sleep(self.delay)
      #Clean up any processes that have now completed
      stillrunning=[]
      for p,fpath in running:
        retcode=p.poll()
        if retcode is None:
          #Process still running
          stillrunning.append((p,fpath))
        else:
          #Process completed
          #Only delete the input file if the process was successful
          if retcode == 0:
            os.remove(str(fpath))
          else:
            print(fpath,retcode)
      #New running list is the list that is still running
      running=stillrunning

#Register for loading from yaml
yaml_manager.register_classes([SimultaneousRequestQueue])

