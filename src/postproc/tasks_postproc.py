
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

def collection(rundict):
  infiles=collect_results.list_result_files(rundict)
  outfpath=osp.join(solnfolder,rundict['master']+'.pkl.gz')
  tdef = {'name':rundict['master'],
          'file_dep':infiles,
          'uptodate':[config_changed(rundict)],
          'targets':[outfpath],
          'actions':[(collect_results.collectall,(rundict,))]}
  return tdef  