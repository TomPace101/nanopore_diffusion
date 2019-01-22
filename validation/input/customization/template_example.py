"""An example of how to provide customization for a template"""

#Function to be bound as a method of the request
def get_template_input(self):
  output={'x': self.x, 'y': y}
  output['z']=self.x+y
  output['n']=self.data['n']
  return output

#List of functions to be bound as methods
request_methods=[get_template_input]

#Function to initialize module variables
def initialize_module(**kwargs):
  for k,v in kwargs.items():
    globals()[k]=v

