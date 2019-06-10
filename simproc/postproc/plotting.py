"""Generate plots"""

#Standard library

#Site packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#This package
from ..requesthandler.customization import CustomizableRequest, make_schema
from ..requesthandler import yaml_manager

FigureRequest_props_schema_yaml="""#FigureRequest
"""

class FigureRequest(CustomizableRequest):
  """
  
  User-defined attributes:
"""
  pass

#Register for loading from yaml
yaml_manager.register_classes([FigureRequest])
