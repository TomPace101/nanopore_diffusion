"""Example calculation function for a job list"""

from simproc.requesthandler.joblist import ExcludeRow

def calc_alpha(self,row,f=0):
  """Calculate """
  alpha = f*(row['b']+row['c']+row['d'])
  if alpha >= 30:
    raise ExcludeRow
  return alpha

request_methods=[calc_alpha]