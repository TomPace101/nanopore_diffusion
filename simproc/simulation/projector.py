"""Project an expression into a function space and save the result."""

#Standard library
from argparse import Namespace

#Site packages
import fenics as fem

#This package
from ..requesthandler import yaml_manager

_ProjectorConditions_props_schema_yaml="""#ProjectorConditions
functionname: {type: string}
functiontype: {type: string}
expression: {type: string}
parameters:
  type: array
  items: {type: string}
projection_kwargs: {type: string}
"""
ProjectorConditions_props_schema=yaml_manager.readstring(_ProjectorConditions_props_schema_yaml)
ProjectorConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    ProjectorConditions_props_schema,['expression'])

class ProjectionSimulator(simrequest.SimulationRequest):
  """Simulator for projecting an expression into a function space

  The result is saved in the ``soln`` attribute.
  
  Conditions:

     - functiontype = optional string specifying type of function:
    
      - 'scalar' (default) for a scalar function
      - 'vector' for a vector function
      
    - expression = required string, the expression to project
    
      Note that the python ``format`` method is called on the string, using the mesh metadata as the keyword arguments.
      This allows the expression to reference variables defining the mesh structure, without using FEniCS parameters.
      
    - parameters = optional list of parameters names to use in the expression, empty for no parameters.
      
      Note that the values of these parameters are taken from the Simulator's ``metadata`` attribute.
      
    - projection_kwargs = optional dictionary of keyworg arguments to the ``project`` function.
    
      This is used, for example, to set the linear solver and preconditioner."""

  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,ProjectorConditions_schema)
    
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Requested Function space
    ftype=getattr(conditions,'functiontype','scalar').lower()
    if ftype == 'scalar':
      self.V = fem.FunctionSpace(self.meshinfo.mesh,'P',self.conditions.elementorder)
    elif ftype == 'vector':
      self.V = fem.VectorFunctionSpace(self.meshinfo.mesh, 'P', self.conditions.elementorder)
    else:
      raise Exception("Invalid functiontype: %s"%ftype)

    #Apply mesh metadata to the expressions
    exprstr=conditions.expression.format(**self.meshinfo.metadata)
    
    #Create the parameters dictionary
    params={}
    for k in getattr(conditions,'parameters',[]):
      params[k]=self.metadata[k]
    
    #Create the expression object
    self.expr=fem.Expression(exprstr,element=self.V.ufl_element(),**params)

    #Get the keyword arguments for projection
    projection_kwargs=getattr(conditions,'projection_kwargs',{})

    #Get the function name
    functionname=getattr(conditions,'functionname','projected')

    #Do the projection
    self.soln=fem.project(self.expr,self.V,**projection_kwargs)

    #Set the function name
    self.soln.rename(functionname,'')

#Register for loading from yaml
yaml_manager.register_classes([ProjectionSimulator])
