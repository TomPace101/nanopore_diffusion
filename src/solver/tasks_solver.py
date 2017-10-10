
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

def dosolve(params):
  assert params.equation=='fickian_unhomog', "Only one equation implemented for now." #TODO: support other equations
  xmlfiles = [x for x in solver_general.List_Mesh_Input_Files(params)] #TODO: this will need to change when we support other equations
  outdir=osp.join(solnfolder,params.meshname)
  outfiles=['conc.pvd','flux.pvd']
  outpaths=[osp.join(outdir,f) for f in outfiles]
  tdef = {'name':params.modelname,
          'file_dep':[osp.join(solverfolder,'fickian_unhomog.py')]+xmlfiles,
          'uptodate':[config_changed(params.__dict__)],
          'targets':outpaths,
          'actions':[(fickian_unhomog.SolveMesh,(params,))]}
  return tdef