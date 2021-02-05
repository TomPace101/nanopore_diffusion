"""Pre/Post-processing calculations for the simulations"""

import fenics as fem

def cell_inverse_integral(self,attrpath,pcell=None,funcpath=None):
  """Compute the inverse of the integral of the inverse
  of the specified function over the given cell.

  Arguments:

    - attrpath = attribute path for storage of result
    - pcell = optional, physical cell number for the cell to calculate integral over,
        empty or None to use entire model
    - funcpath = optional, path to function to integrate,
        empty or None to use a constant value of 1.0

  Required attributes:

    - meshinfo.mesh = FEniCS Mesh object
    - meshinfo.cells = FEniCS MeshFunction object for cell numbers
  
  No return value."""
  if funcpath is None:
    this_func=fem.Constant(1.0)
  else:
    this_func=self.get_nested(funcpath)
  if pcell is None:
    this_dx=fem.Measure('cell',domain=self.meshinfo.mesh)
  else:
    this_measure=fem.Measure('cell',domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.cells)
    this_dx=this_measure(pcell)
  result=1.0/fem.assemble((1.0/this_func)*this_dx)
  self.set_nested(attrpath,result)

#List of functions to be bound as methods
request_methods=[cell_inverse_integral]
