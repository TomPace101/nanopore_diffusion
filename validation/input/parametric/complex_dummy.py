"""An example of customization for request generation"""


import sys
sys.path.append('../../..')

import simproc.requesthandler.locators as locators

def get_child_kwargs(self,index,prefix,x,y):
  """Compute a keyword arguments dictionary from the input dictionary"""
  out={}
  chnum='%03d'%index
  out['name']='parametric.customized.%s'%chnum
  out['test']='%s %s.\n'%(prefix,chnum)
  out['test']+="x=%d, y=%d"%(x,y)
  out['outfile']=locators.OutputFile("child_%s.txt"%chnum)
  return out

#List of functions to be bound as methods
add_methods=[get_child_kwargs]
