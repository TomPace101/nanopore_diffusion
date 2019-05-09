"""Solve the multi-species unhomogenized Smoluchowski diffusion problem,
with a reactive boundary condition,
and extract the data needed for post-processing efforts.

Species may vary within each subdomain."""

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

_SUConditions_props_schema_yaml="""#SUConditions
species: {type: object}
"""
SUConditions_props_schema=yaml_manager.readstring(_SUConditions_props_schema_yaml)
SUConditions_schema=simrequest.update_schema_props(simrequest.GenericConditions_schema,
                                                    SUConditions_props_schema,['species'])

class SUSimulator(simrequest.SimulationRequest):
  """Simulator for for Unhomogenized Smoluchowski Diffusion"""
  
  _props_schema=simrequest.update_conditions(simrequest.SimulationRequest._props_schema,SUConditions_schema)
    
  def run_sim(self):

    #For convenience
    conditions=Namespace(**self.conditions)

#Register for loading from yaml
yaml_manager.register_classes([SUSimulator])

