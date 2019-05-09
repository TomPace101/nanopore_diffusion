"""Solve the linearized Poisson-Boltzmann equation"""

#Standard library
from argparse import Namespace

#Site packages
import numpy as np
import fenics as fem

#This package
from ..requesthandler import yaml_manager
from .meshinfo import MeshInfo
from . import simrequest
from . import equationbuilder

_LPBConditions_props_schema_yaml="""#LPBConditions
kappa: {type: number}
"""
LPBConditions_props_schema=yaml_manager.readstring(_LPBConditions_props_schema_yaml)
LPBConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    LPBConditions_props_schema,['kappa'])

class LPBSimulator(simrequest.SimulationRequest):
  """Simulator for linearized Poisson-Boltzmann equation
  
  User-defined attributes:
  
    - """
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,LPBConditions_schema)    
  
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Properties of problem domain
    self.kappa = self.conditions.kappa
    
    


#Register for loading from yaml
yaml_manager.register_classes([LPBSimulator])
