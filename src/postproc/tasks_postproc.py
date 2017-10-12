
#Standard library
import os
import os.path as osp
import sys

#Site packages
from doit.tools import config_changed

#Local
from folderstructure import *
import useful
import collect_results

def collection(run):
  infiles=collect_results.list_result_files(run)
  outfpath=osp.join(solnfolder,rundict['master']+'.pkl.gz')
  tdef = {'name':run.master,
          'file_dep':infiles,
          'uptodate':[config_changed(run.to_dict())],
          'targets':[outfpath],
          'actions':[(collect_results.collectall,(run,))]}
  return tdef  