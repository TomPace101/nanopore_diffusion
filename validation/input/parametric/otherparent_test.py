"""An example of using another parent request to define a parametric request"""

import sys
sys.path.append('../../..')

import simproc.requesthandler.locators as locators

def get_child_kwargs(self,index,prefix,z,xy_request):
  """Compute a keyword arguments dictionary from the input dictionary"""
  out={}
  chnum='%03d'%index
  out['name']='parametric.with_other_parent.%s'%chnum
  out['test']='%s %s.\n'%(prefix,chnum)
  out['test']+="z=%d\n"%(z)
  out['test']+="The other request contained test data:\n"
  for line in xy_request.test.split('\n'):
    out['test']+="\t%s\n"%line
  out['outfile']=locators.OutputFile("zchild_%s.txt"%chnum)
  return out

#List of functions to be bound as methods
request_methods=[get_child_kwargs]
