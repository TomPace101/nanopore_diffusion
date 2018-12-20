"""An example of using another parent request to define a parametric request"""

import sys
sys.path.append('../../..')

import simproc.requesthandler.locators as locators

def get_child_kwargs(self,fields):
  """Compute a keyword arguments dictionary from the input dictionary"""
  out={}
  chnum='%03d'%fields['index']
  out['name']='parametric.customized.%s'%chnum
  out['test']=fields['prefix']+' '+chnum+'.\n'
  out['test']+="z=%d\n"%(fields['z'])
  out['test']+="The other request contained test data:\n"
  for line in fields['xy_request'].test.split('\n'):
    out['test']+="\t%s\n"%line
  out['outfile']=locators.OutputFile("zchild_%s.txt"%chnum)
  return out

#List of functions to be bound as methods
add_methods=[get_child_kwargs]
