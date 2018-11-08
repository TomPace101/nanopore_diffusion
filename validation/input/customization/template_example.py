"""An example of how to provide customization for a template"""

def get_template_input(self):
  output={'x': self.x, 'y': y}
  output['z']=self.x+y
  return output

def initialize_module(**kwargs):
  for k,v in kwargs.items():
    globals()[k]=v
