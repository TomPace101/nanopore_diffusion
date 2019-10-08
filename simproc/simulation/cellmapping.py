"""Project an expression into a function space and save the result."""

#Standard library
from argparse import Namespace

#Site packages
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from . import simrequest
from . import expressions


_MappingConditions_props_schema_yaml="""#MappingConditions
functionname: {type: string}
functiontype: {type: string}
mapping: {type: object}
projection_kwargs: {type: string}
"""
MappingConditions_props_schema=yaml_manager.readstring(_MappingConditions_props_schema_yaml)
MappingConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    MappingConditions_props_schema,['mapping'])

class CellMappingSimulator(simrequest.SimulationRequest):
  """Simulator for mapping mesh function cell values to a function in the given function space

  The result is saved in the ``soln`` attribute.
  
  Conditions:

    - functionname = string specifying name for projected function (not the attribute)

    - functiontype = optional string specifying type of function:
    
      - 'scalar' (default) for a scalar function
      - 'vector' for a vector function
      - 'matrix' for a rank-2 tensor function
      
    - mapping = required dictionary, the mapping from cell value to function value
    
    - projection_kwargs = optional dictionary of keyword arguments to the ``project`` function.
    
      This is used, for example, to set the linear solver and preconditioner."""

  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,MappingConditions_schema)
    
  def run_sim(self):

     #For convenience
    conditions=Namespace(**self.conditions)
    spatial_dims=self.meshinfo.mesh.geometry().dim()

    #Requested Function space, and the UserExpression subclass instance
    ftype=getattr(conditions,'functiontype','scalar').lower()
    if ftype == 'scalar':
      self.V = fem.FunctionSpace(self.meshinfo.mesh,'P', conditions.elementorder)
    elif ftype == 'vector':
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder)
    elif ftype == 'matrix':
      self.V = fem.TensorFunctionSpace(self.meshinfo.mesh, 'P', conditions.elementorder, (spatial_dims,spatial_dims))
    else:
      raise Exception("Invalid functiontype: %s"%ftype)

    #Create the expression
    self.loadcellmapping('expr',conditions.mapping,ftype)

    #Get the keyword arguments for projection
    projection_kwargs=getattr(conditions,'projection_kwargs',{})

    #Get the function name
    functionname=getattr(conditions,'functionname','projected')

    #Do the projection
    self.soln=fem.project(self.expr,self.V,**projection_kwargs)

    #Set the function name
    self.soln.rename(functionname,'')

#Register for loading from yaml
yaml_manager.register_classes([CellMappingSimulator])
