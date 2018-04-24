#Standard library
from __future__ import print_function, division #Python 2 compatibility
import importlib
import os.path as osp
import sys

#Site packages

#Local
import folderstructure as FS
import common
import buildgeom
import solver_general

#Solver modules
solver_module_list=['fickian_unhomog','smol_unhomog','tdpnp_unhomog']

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

class ModelParameters(solver_general.ModelParametersBase):
  """Extend ModelParametersBase to allow for task generation and execution of arbitrary solver module
  Attributes:
    solverclass = the class used to solve the model"""
  __slots__=('solverclass',)
  #Load the solver modules and map equation names to their solver classes
  #Mapping from ModelParameters.equation to the appropriate solver classes
  #Each module implementing a solver should define a solverclasses dictionary, so we just need to put them all together
  solverclasses={}
  for sm_name in solver_module_list:
    sm=importlib.import_module(sm_name)
    solverclasses.update(sm.solverclasses)
  
  def __init__(self,**kwd):
    #Initialization from base class
    super(ModelParameters, self).__init__(**kwd)
    #Get the solver class
    self.solverclass=self.solverclasses[self.equation]
    #Add to code file dependencies
    self._more_inputfiles+=[thisfile, sys.modules[self.solverclass.__module__].__file__]

  def run(self):
    """Run the loaded model."""
    #Run the solver
    print(self.modelname)
    self.solverclass.complete(self)

#Support command-line arguments
if __name__ == '__main__':
  program_description='Solve a diffusion equation with fenics'
  input_file_description='Path to file containing ModelParameters definitions'
  other_selection={'equation':ModelParameters.solverclasses.keys()}
  
  common.run_cmd_line(program_description,input_file_description,ModelParameters,other_selection=other_selection)
  #other_selection is needed so we only try to run models whose equation we have a solver for
