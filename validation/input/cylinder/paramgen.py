"""Generate the variations of the cylinder analysis"""

def get_template_input(self):
  output={}
  paramslist=[('%03d'%idx,r) for idx,r in enumerate(self.data['r_values'])]
  output['paramslist']=paramslist
  return output

#List of functions to be bound as methods
request_methods=[get_template_input]
