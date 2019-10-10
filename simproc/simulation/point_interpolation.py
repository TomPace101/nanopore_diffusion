"""Create a function by interpolating from point data in a CSV file"""

#Standard library
from argparse import Namespace

#Site packages
import numpy as np
import fenics as fem
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

#This package
from ..requesthandler import yaml_manager
from ..requesthandler import timing
from . import simrequest
from . import meshinfo

_InterpolationConditions_props_schema_yaml="""#InterpolationConditions
functionname: {type: string}
boundaryvalue: {type: number}
"""
InterpolationConditions_props_schema=yaml_manager.readstring(_InterpolationConditions_props_schema_yaml)
InterpolationConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    InterpolationConditions_props_schema,[])

class InterpolationSimulator(simrequest.SimulationRequest):
  """Simulator for projecting an expression into a function space

  The point data must be loaded with a command in ``loaddata``.
  The result of the interpolation is saved in the ``soln`` attribute.
  
  Conditions:

    - functionname = string specifying name for projected function (not the attribute)
    - boundaryvalue = function value to be used at mesh boundaries"""

  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,InterpolationConditions_schema)
    
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)
    d=self.meshinfo.spatial_dims()
    meshcoords=self.meshinfo.coordinates()
    Nmesh=meshcoords.shape[0]

    #Function space, and its coordinates
    self.V = fem.FunctionSpace(self.meshinfo.mesh,'P', conditions.elementorder)
    Ndof_pts=self.V.dim()
    dofinfo=meshinfo.DOFInfo(self.meshinfo,self.V)
    dofcoords=dofinfo.coordinates()

    #Load the input data
    self.process_load_commands()

    #Get the function name
    functionname=getattr(conditions,'functionname','projected')

    #Get boundary and non-boundary DOF coordinates
    boundary_pts, nonbound_pts = dofinfo.boundary_points()

    #Array of boundary values
    boundary_vals=np.array([conditions.boundaryvalue]*boundary_pts.shape[0])



  def loadcsv(self,infpath,attrpath="pointdata"):
    """Load a CSV file containing scalar function values at different points

    Arguments:

      - infpath = path to the input CSV file
      - attrpath = attribute path to load the data into
    
    No return value."""
    fpt=self.renderstr(self.get_stored(infpath))
    df=pd.read_csv(fpt)
    self.set_nested(attrpath,df)
    return

#Register for loading from yaml
yaml_manager.register_classes([InterpolationSimulator])
