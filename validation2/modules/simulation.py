"""Pre/Post-processing calculations for the simulations"""

import numpy as np
from scipy.integrate import trapz
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

def calc_ratio(self, numerator, denominator, outattr, use_stored=True):
  if use_stored:
    a=self.get_stored(numerator)
    b=self.get_stored(denominator)
  else:
    a=self.get_nested(numerator)
    b=self.get_nested(denominator)
  res=a/b
  self.set_nested(outattr, res)

def calc_product(self,factor1,factor2,outattr):
  a=self.get_stored(factor1)
  b=self.get_stored(factor2)
  res=a*b
  self.set_nested(outattr,res)

def calc_delta(self, vmin, vmax, outattr):
  """Compute a delta from min and max

  Arguments:
  
    - vmin = minimum value
    - vmax = maximum value
    - outattr = attribute path for storing result"""
  dv=self.get_stored(vmax)-self.get_stored(vmin)
  self.set_nested(outattr,dv)
  return

def project_exp_pot(self,outattr="exp_pot",funcname="exp_pot",solver="cg",precond="amg",):
  """Project the exponential of the potential, for the calculation of Xi

  Arguments:

    - outattr: attribute path for storing the result
    - funcname: name to be given to the projected function
    - solver: linear solver to be used for the projection operation
    - precond: preconditioner to be used for the projection operation

  """
  beta = self.conditions['beta']
  expr = fem.exp(-beta*self.potential)
  res = fem.project(expr,self.scalar_V,solver_type=solver,preconditioner_type=precond)
  res.rename(funcname,funcname)
  self.set_nested(outattr,res)

def calc_expected_result(self,outattr):
  """This is only for the PMF x-variation case"""
  xvals=np.linspace(0.0,0.5,1000,endpoint=True)
  yvals=np.exp(xvals**2)
  intg=2*trapz(yvals,xvals)
  res=float(np.exp(0.25)/intg)
  self.set_nested(outattr,res)

#List of functions to be bound as methods
request_methods=[cell_inverse_integral, calc_ratio, calc_product, calc_delta, project_exp_pot, calc_expected_result]
