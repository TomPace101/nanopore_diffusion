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
coordcolumns:
  type: array
  items: {type: string}
valuecolumn: {type: string}
boundaryvalue: {type: number}
"""
InterpolationConditions_props_schema=yaml_manager.readstring(_InterpolationConditions_props_schema_yaml)
InterpolationConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    InterpolationConditions_props_schema,
                                                    ['boundaryvalue'])

class InterpolationSimulator(simrequest.SimulationRequest):
  """Simulator for projecting an expression into a function space

  The point data must be loaded with a command in ``loaddata``,
  and saved in the attribute ``pointdata``.
  The result of the interpolation is saved in the ``soln`` attribute.
  
  Conditions:

    - functionname = optional string specifying name for projected function (not the attribute)
    - boundaryvalue = required function value to be used at mesh boundaries
    - coordcolumns = optional list of names (as strings) for the columns containing X,Y,and Z, respectively.
    
      If there are only two dimensions, just leave off the column for Z.
      The default value is [x,y,z].

    - valuecolumn = optional name, as string, of the column containing the function value. (Defaults to 'f')

      Note that column names are all case-sensitive."""

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
    df=self.pointdata

    #Get the function name
    functionname=getattr(conditions,'functionname','projected')

    #Get the column names for the input data
    coordcolumns=getattr(conditions,'coordcolumns',['x','y','z'])
    valuecolumn=getattr(conditions,'valuecolumn','f')

    #Get boundary and non-boundary DOF coordinates
    boundary_pts, nonbound_pts = dofinfo.boundary_points()

    #Array of boundary values
    boundary_vals=np.array([conditions.boundaryvalue]*boundary_pts.shape[0])

    #Input coordinates and function values, as separate arrays
    inpts=df.loc[:,coordcolumns].values
    invals=df.loc[:,[valuecolumn]].values.flatten()

    #Combine boundary and input points
    data_pts=np.vstack([inpts,boundary_pts])
    data_vals=np.hstack([invals,boundary_vals])

    #Set up the interpolator
    self.interp_setup_timer=timing.Timer()
    ilator=LinearNDInterpolator(data_pts,data_vals,fill_value=conditions.boundaryvalue,rescale=False)
    self.interp_setup_timer.stop()
    #Interpolate at all dof points
    self.interp_run_timer=timing.Timer()
    dofvals=ilator(dofcoords)
    self.interp_run_timer.stop()

    #Define function from the interpolated dof values
    self.soln=fem.Function(self.V,name=functionname)
    junk=self.soln.vector().set_local(dofvals)
    del junk

    #Done
    return

  def compute_residual_errors(self,dfpath="pointdata",funcattr="soln",outattr="residuals"):
    df=self.get_nested(dfpath)
    func=self.get_nested(funcattr)
    #Get the column names for the input data
    coordcolumns=self.get_nested_default('conditions.coordcolumns',['x','y','z'])
    valuecolumn=self.get_nested_default('conditions.valuecolumn','f')
    #Set up the output dataframe
    resid=df.copy()
    #Compute the function value at each input point
    def calcvalue(row):
      pt=[row[c] for c in coordcolumns]
      return func(*pt)
    resid['interpolated']=df.apply(calcvalue,axis=1)
    #Compute the residuals
    def calcresidual(row):
      return row['interpolated']-row[valuecolumn]
    resid['error']=resid.apply(calcresidual,axis=1)
    #Store result
    self.set_nested(outattr,resid)
    return

#Register for loading from yaml
yaml_manager.register_classes([InterpolationSimulator])
