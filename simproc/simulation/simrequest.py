"""Base functionality for simulation requests"""

#Site packages
import fenics as fem

#This package
from ..requesthandler.customization import CustomizableRequest
from .meshinfo import MeshInfo

_SimulationRequest_props_schema_yaml="""#SimulationRequest
mesh:
  anyOf:
    - type: pathlike
    - type: object
conditions:
  type: object
dataextraction:
  type: array
loaddata:
  type: array
metadata:
  type: object
"""

class SimulationRequest(CustomizableRequest):
  """Base class for FEniCS simulations
  
  User-defined attributes:
  
    - mesh: path to mesh hdf5 file, or dictionary specifying mesh, depending on the simulation
    - conditions: dictionary specifying model conditions such as element order, boundary conditions, etc.
    - dataextraction: TODO
    - loaddata: TODO
    - metadata: TODO"""##TODO
  _self_task=True
  _required_attrs=[] ##TODO
  _outputfile_attrs=[] ##TODO
  _inputfile_attrs=[] ##TODO
  _props_schema=readyaml(_SimulationRequest_props_schema_yaml)

#Register for loading from yaml
register_classes([SimulationRequest])

