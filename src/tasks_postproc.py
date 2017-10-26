
#Standard library
import os
import os.path as osp

#Site packages

#Local
from folderstructure import *
import useful
import collect_results

def collection(basename,model_list):
  infiles=collect_results.list_inputfiles(basename,model_list)
  outfpath=osp.join(solnfolder,basename+'.pkl.gz')
  tdef = {'name': basename,
          'file_dep':infiles+[collect_results.__file__],
          'targets':[outfpath],
          'actions':[(collect_results.do_collection,(infiles,outfpath))]}
  return tdef  