
#Standard library
import os
import os.path as osp
import sys

#Site packages
from doit.tools import config_changed

#Local
sys.path.append('..')
import useful
from folderstructure import *
from . import fickian_unhomog

def dosolve(params):
  xmlfiles = [x for x in fickian_unhomog.List_Mesh_Input_Files(params)]
  outdir=osp.join(solndir,params.meshname)
  outfiles=['conc.pvd','flux.pvd']
  outpaths=[osp.join(outdir,f) for f in outfiles]
  tdef = {'name':params.bcname,
          'file_dep':[osp.join(solverfolder,'fickian_unhomog.py')]+xmlfiles,
          'uptodate':[config_changed(params.__dict__)],
          'targets':outpaths,
          'actions':[(fickian_unhomog.SolveMesh,(params,))]}
  return tdef