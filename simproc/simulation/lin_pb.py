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

_LPBSimulator_props_schema_yaml="""#LPBSimulator
"""

class LPBSimulator(simrequest.SimulationRequest):
  """Simulator for linearized Poisson-Boltzmann equation
  
  User-defined attributes:
  
    - """
    
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

    #Properties of problem domain
    self.kappa = self.conditions.kappa
    
    


#Register for loading from yaml
yaml_manager.register_classes([LPBSimulator])
