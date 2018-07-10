"""Solve the homogenized Fickian diffusion problem
and extract the data needed for post-processing efforts"""

#Standard library

#Site packages
import fenics as fem

#Local
import unitsystem as UN
import common
import simulator_general

class EquationTerm(simulator_general.EquationTermBase):
  __slots__=('bilinear')

class Conditions(simulator_general.GenericConditions):
  """Condition defnitions for use with HomogFickianSimulator
  
  Additional attributes:
  
    - boundaries = list of physical lines representing boundaries, for construction of boundary terms"""
  __slots__=('boundaries')

class PeriodicBoundary(fem.SubDomain):
  """SubDomain subclass for Periodic boundary condition"""
  def __init__(self,xlims,ylims):
    """Arguments:
    
      - xlims = pair of x-values: (xmin,xmax)
      - ylims = pair of y-values: (ymin,ymax)"""
    super(PeriodicBoundary, self).__init__()
    self.left,  self.right = xlims
    self.bottom,self.top   = ylims
    self.xspan = self.right-self.left
    self.yspan = self.top-self.bottom
  # Left boundary is "target domain" G
  def inside(self, x, on_boundary):
    # return True if on left or bottom boundary AND NOT on one of the two corners 
    #  (self.left, self.top) and (self.right, self.bottom)
    return bool((fem.near(x[0], self.left) or fem.near(x[1], self.bottom)) and 
                (not ((fem.near(x[0], self.left) and fem.near(x[1], self.top)) or 
                      (fem.near(x[0], self.right) and fem.near(x[1], self.bottom)))) and on_boundary)
  def map(self, x, y):
    if fem.near(x[0], self.right) and fem.near(x[1], self.top):
      y[0] = x[0] - self.xspan
      y[1] = x[1] - self.yspan
    elif fem.near(x[0], self.right):
      y[0] = x[0] - self.xspan
      y[1] = x[1]
    else:   # fem.near(x[1], self.top)
      y[0] = x[0]
      y[1] = x[1] - self.yspan

class HomogFickianSimulator(simulator_general.GenericSimulator):
  """Simulator for Homogenized Fickian Diffusion

  Additional attributes not inherited from GenericSimulator:

    - meshinfo = instance of simulator_general.MeshInfo
    - conditions = instance of HomogFickianConditions
"""
  def __init__(self,modelparams):
    """Initialize the model.

    Arguments:

      - modelparams = simulator_run.ModelParameters instance"""

    #Load parameters, init output, load mesh
    super(HomogFickianSimulator, self).__init__(modelparams)


simulatorclasses={'fickian_homog':HomogFickianSimulatorSimulator}