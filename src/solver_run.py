#Standard library
import importlib

#Site packages

#Local
import useful
import solver_general

#Solver modules
solver_module_list=['fickian_unhomog','smol_unhomog']

#Load solver classes from the solver modules
#Mapping from solver_general.ModelParameters.equation to the appropriate solver classes
#Each module implementing a solver should define a solverclasses dictionary, so we just need to put them all together
solverclasses={}
for sm_name in solver_module_list:
  sm=importlib.import_module(sm_name)
  solverclasses.update(sm.solverclasses)

#Support command-line arguments
if __name__ == '__main__':
  program_description='Solve a diffusion equation with fenics'
  input_file_description='Path to file containing ModelParameters definitions'
  other_selection={'equation':solverclasses.keys()}
  
  useful.run_cmd_line(program_description,input_file_description,
    solver_general.ModelParameters,
    solver_general.complete_by_ModelParameters,
    other_selection,
    [solverclasses])
