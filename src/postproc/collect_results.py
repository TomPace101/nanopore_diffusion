
#Standard library
import os
import os.path as osp

#Site packages
import pandas as pd

#Local
from folderstructure import *
import useful

def read_all_dicts(rundict):
  #Dictionary of all meshes, by name
  allmeshes={}
  for doc in useful.readyaml_multidoc(rundict['meshparams']):
    allmeshes[doc['meshname']]=doc
  #Dictionary of all models, by name
  allmodels={}
  for doc in useful.readyaml_multidoc(rundict['modelparams']):
    allmodels[doc['modelname']]=doc
  return allmeshes, allmodels

def list_result_files(rundict):
  allmeshes, allmodels=read_all_dicts(rundict)
  infiles=[osp.join(solnfolder,mname,'results.yaml') for mname in allmodels.keys()]
  return infiles

def gen_dataframe(rundict):
  infiles=list_result_files(rundict)
  alldicts=[]
  for infpath in infiles:
    alldicts.append(useful.readyaml(infpath))
  df = pd.DataFrame(alldicts)
  return df


def collectall(rundict):
  outfpath=osp.join(solnfolder,rundict['master']+'.pkl.gz')
  df = gen_dataframe(rundict)
  df.to_pickle(outfpath)