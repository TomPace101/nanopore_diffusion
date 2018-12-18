"""An example of customization for request generation"""


import sys
sys.path.append('../../..')

import simproc.requesthandler.locators as locators

def get_child_kwargs(self,fields):
  """Compute a keyword arguments dictionary from the input dictionary"""
  out={}
  chnum='%03d'%fields['index']
  out['name']='parametric.customized.%s'%chnum
  out['test']=fields['prefix']+' '+chnum+'.\n'
  out['test']+="x=%d, y=%d"%(fields['x'],fields['y'])
  out['outfile']=locators.OutputFile("child_%s.txt"%chnum)
  return out

#List of functions to be bound as methods
add_methods=[get_child_kwargs]
