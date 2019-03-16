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
tmpfile: {type: pathlike}
"""

class MPIRunRequest(request.Request):
  """Run the child request using mpirun
  
  User-defined attributes:
  
    - numproc: number of processes to specify with mpirun
    - child: a request to run
    - tmpfile: optional, Path to temporary file to contain the request"""
  _self_task=False
  _required_attrs=['numproc','child']
  _props_schema=request.make_schema(_MPIRunRequest_props_schema_yaml)
  _child_attrs=['child']
  def __init__(self,**kwargs):
    #Initialization from base class
    super(MPIRunRequest, self).__init__(**kwargs)
    if hasattr(self,'tmpfile'):
      self.tmpfile=filepath.Path(self.tmpfile).expanduser().resolve()
    else:
      self.tmpfile=filepath.Path('.').expanduser().resolve() / 'tmp.yaml'
  def run(self):
    #Create the output directories for the request, one time, to avoid conflicts
    self.child.assure_output_dirs()
    #Write the input file
    yaml_manager.writefile([self.child],str(self.tmpfile))
    #Call MPIrun to start the processes
    args=('mpirun','-np','%d'%self.numproc,sys.executable,'-m',main_mod_name,str(self.tmpfile))
    p=subprocess.Popen(args,cwd=work_path,shell=False)
    #Wait for completion
    retcode=p.wait()
    #Clean up
    os.remove(str(self.tmpfile))
    
#Register for loading from yaml
yaml_manager.register_classes([MPIRunRequest])
