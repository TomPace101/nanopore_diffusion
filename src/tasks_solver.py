
#Standard library
import os
import os.path as osp
import sys

#Site packages
from doit.tools import config_changed

#Local
from folderstructure import *
import solver_general
import fickian_unhomog

#Mapping from solver_general.ModelParameters.equation to the appropriate solver classes
solverclasses={'fickian_unhomog':fickian_unhomog.UnhomogFickianSolver,
               'smol_unhomog':smol_unhomog.SUSolver}

def initobj(c,*args):
  "A doit action can't just initialize an object, so this function does the init"
  obj=c(*args)
  return

def list_outputfiles(cmdlist):
  """Get a list of all the files generated by the data extraction commands.
  The info.yaml file is included as well.
  Arguments:
    cmdlist = list of data extraction commands,
      each command consists of pair (cmdname, arguments):
        cmdname = name of data extraction method of the solver class
        arguments = dictionary of all arguments needed by the extraction method
  Return:
    outfiles = list of generated output files (names only, not including their folder)"""
  #Currently, we assume all files can only come from the 'filename' argument
  filearg_list=['filename']
  outfiles=['info.yaml']
  for cmdname, arguments in cmdlist:
    #Check all possible arguments that could contain the name of an output file
    present_args=[n for n in filearg_list if n in arguments.keys()]
    outfiles.extend([arguments[n] for n in present_args])
  return outfiles

def dosolve(modelparams,meshparams):
  assert modelparams.equation in solverclasses.keys(), "Unrecognized equation: %s"%modelparams.equation
  solver_class=solverclasses[modelparams.equation]
  solver_codefile=sys.modules[solver_class.__module__].__file__
  xmlfiles = [x for x in solver_general.List_Mesh_Input_Files(modelparams.meshname,meshparams.basename)]
  outdir=osp.join(solnfolder,modelparams.basename,modelparams.modelname)
  outfiles=list_outputfiles(modelparams.dataextraction)
  outpaths=[osp.join(outdir,f) for f in outfiles]
  tdef = {'name':modelparams.modelname,
          'file_dep':[solver_general.__file__,osp.join(srcfolder,solver_codefile)]+xmlfiles, #This takes care of dependency on mesh
          'uptodate':[config_changed(modelparams.to_dict())], #Mesh dependency already taken care of in file_dep
          'targets':outpaths,
          'actions':[(initobj,(solver_class,modelparams,meshparams,True))]}
  return tdef