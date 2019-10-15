"""Generate the variations of the cylinder analysis"""

def get_template_input(self):
  output={}
  output['r_list']=[('%02d'%idx,r) for idx,r in enumerate(self.data['r_values'])]
  output['pot_list']=[('%01d'%idx,p) for idx,p in enumerate(self.data['pot_values'])]
  output['collection_files']=[]
  return output

#List of functions to be bound as methods
request_methods=[get_template_input]
