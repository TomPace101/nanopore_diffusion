"""Run requests with mpirun"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import subprocess
import sys
import os

#This package
from . import filepath
from . import request
from . import yaml_manager
from . import simultaneous

#Constants from simultaneous
main_mod_name=simultaneous.main_mod_name
work_path=simultaneous.work_path

_MPIRunRequest_props_schema_yaml="""#MPIRunRequest
name: {type: string}
numproc:
  type: integer
  minimum: 1
child: {type: request}
outfile: {type: pathlike}
errfile:
  anyOf:
    - {type: pathlike}
    - {type: 'null'}
tmpfile: {type: pathlike}
"""

class MPIRunRequest(request.Request):
  """Run the child request using mpirun
  
  User-defined attributes:
  
    - numproc: number of processes to specify with mpirun
    - child: a request to run
    - outfile: Path to output file
    - errfile: optional, Path to error output file, or None to redirect to `outfile`
    - tmpfile: optional, Path to temporary file to contain the request"""
  _self_task=False
  _required_attrs=['numproc','child','outfile']
  _props_schema=request.make_schema(_MPIRunRequest_props_schema_yaml)
  _child_attrs=['child']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(SimultaneousRequestQueue, self).__init__(**kwargs)
    if hasattr(self,'tmpfile'):
      self.tmpfile=filepath.Path(self.tmpfile).expanduser().resolve()
    else:
      self.tmpfile=filepath.Path('.').expanduser().resolve() / 'tmp.yaml'
    if hasattr(self,'errfile'):
      self.errfile=filepath.Path(self.errfile).expanduser().resolve()
    else:
      self.errfile=None
  def run(self):
    #Write the input file
    yaml_manager.writefile([req],str(self.tmpfile))
    #Call MPIrun to start the processes
    args=('mpirun','-np','%d'%,sys.executable,'-m',main_mod_name,str(self.tmpfile))
    stderr = subprocess.STDOUT if self.errfile is None else str(self.errfile)
    p=subprocess.Popen(args,cwd=work_path,shell=False,stdout=str(self.outfile),stderr=stderr)
    #Wait for completion
    retcode=p.wait()
    #Clean up
    os.remove(str(self.tmpfile))
    
#Register for loading from yaml
yaml_manager.register_classes([MPIRunRequest])
