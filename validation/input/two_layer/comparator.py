"Compare simulation result to analytical result"

from simproc.requesthandler.nested import Stored

def calc_exact(self,attrpath="expected_answer"):
  """Calculate the theoretically correct answer"""
  #Get input values
  alpha1=self.meshinfo.metadata['alpha1']
  alpha2=self.meshinfo.metadata['alpha2']
  D1=self.metadata['D1']
  D2=self.metadata['D2']
  #Calculations
  parval=D1*alpha1 + D2*alpha2
  perpval=D1*D2/(alpha1*D2+alpha2*D1)
  ans=[[parval,0],[0,perpval]]
  #Store result
  self.set_nested(attrpath,ans)
  return

def compare_answers(self,summarypath="comparison_ok",diffpath="error_values",expected="expected_answer",received="received_answer",difftol=0.01):
  #Get input values
  ans_correct=self.get_stored(expected)
  ans_found=self.get_stored(received)
  tolerance=self.get_stored(difftol)
  #Compute the differences and compare to the tolerance
  all_ok=True
  err_vals=[]
  for i,corr_row in enumerate(ans_correct):
    err_row=[]
    for j,corr_value in enumerate(corr_row):
      diff=corr_value-ans_found[i][j]
      err_row.append(diff)
      all_ok = all_ok and abs(diff)<=tolerance
    err_vals.append(err_row)
  #Store result
  self.set_nested(summarypath,all_ok)
  self.set_nested(diffpath,err_vals)
  return


#List of functions to be bound as methods
request_methods=[calc_exact, compare_answers]
