
#Standard library
import os
import os.path as osp
import sys

#Site packages
from doit.tools import config_changed

#Local
from folderstructure import *
import useful
import solver_general
import fickian_unhomog

solverfuncs={'fickian_unhomog':fickian_unhomog.SolveMesh}

def dosolve(modelparams,meshparams):
  assert modelparams.equation in solverfuncs, "Unrecognized equation: %s"%modelparams.equation
  xmlfiles = [x for x in solver_general.List_Mesh_Input_Files(modelparams.meshname)]
  outdir=osp.join(solnfolder,modelparams.modelname)
  outfiles=['results.yaml'] #TODO: there can be other files output as well, depending on what data is extracted
  outpaths=[osp.join(outdir,f) for f in outfiles]
  tdef = {'name':modelparams.modelname,
          'file_dep':[osp.join(solverfolder,'fickian_unhomog.py')]+xmlfiles,
          'uptodate':[config_changed(params.to_dict())],
          'targets':outpaths,
          'actions':[(fickian_unhomog.SolveMesh,(modelparams,meshparams))]}
  return tdef