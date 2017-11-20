
#Standard library
import os
import os.path as osp

#Site packages
import pandas as pd

#Local
from folderstructure import *
import useful

#Constants
infofile='info.yaml'
collected_df_fname='collected_results.pkl.gz'

def list_inputfiles_foldersearch(basename):
  """Return a list of all info.yaml files for the given basename.
  The list is generated by reading the directories.
  This won't work if the analyses haven't been done yet."""
  basedir=osp.join(solnfolder,basename)
  assert osp.isdir(basedir), "Could not find base directory: %s"%basedir
  filelist=[]
  for content in os.listdir(basedir):
    dirname = osp.join(basedir,content)
    if osp.isdir(dirname):
      fpath=osp.join(dirname,infofile)
      if osp.isfile(fpath):
        filelist.append(fpath)
  return filelist

def list_inputfiles(basename,model_list):
  """Return a list of all info.yaml files for the given basename.
  The list is generated by looking at all the models."""
  basedir=osp.join(solnfolder,basename)
  filelist=[]
  for modelparams in model_list:
    if modelparams.basename==basename:
      #Add this one to the list
      filelist.append(osp.join(basedir,modelparams.modelname,'info.yaml'))
  return filelist

def get_columns(d,exclusions):
  """Return list of all keys where the values are integers, floats, or strings,
  and call recursively on any value that has its own 'items' attribute."""
  cols=[]
  for k,v in d.items():
    if not k in exclusions:
      if type(v)==int or type(v)==float or type(v)==str:
        if not k in cols:
          cols.append(k)
      elif hasattr(v,'items'):
        newcols=get_columns(v,exclusions)
        for c in newcols:
          if c not in cols:
            cols.append(c)
  return cols

def get_all_columns(dlist,exclusions):
  """Return the superset of columns for each d in dlist"""
  columns=get_columns(dlist[0],exclusions)
  for d in dlist[1:]:
    newcols=get_columns(d,exclusions)
    for c in newcols:
      if not c in columns:
        columns.append(c)
  return columns

def flatdict(d,cols,name,exclusions):
  """Flatten a potentially nested dictionary so it can be added to a DataFrame with the specified columns."""
  fd={}
  for k,v in d.items():
    if not k in exclusions:
      if k in cols and (type(v)==int or type(v)==float or type(v)==str):
        fd[k]=v
      elif hasattr(v,'items'):
        newname=name+'->'+k
        for newk, newv in flatdict(v,cols,newname,exclusions).items():
          if newk in fd.keys():
            assert newv==fd[newk], "In %s, unequal assignments to %s: previously had %s, but %s wants to change it to %s"%(name,str(newk),str(fd[newk]),str(k),str(newv))
          else:
            fd[newk]=newv
  return fd

def dicts_to_dataframe(alldicts,exclusions):
  """Create a pandas dataframe from an iterable of dictionaries.
  Arguments:
    alldicts = iterable of dictionaries
    exclusions = list of keys to exlcude
  Return value:
    df = pandas dataframe
  For each dictionary:
    Anything that is a number or string is added directly.
    Anything that is a dictionary has its items treated the same way.
      (More specifically, anything that has an 'items' attribute.)
    Everything else is ignored, including any sequences."""
  #Get the list of columns for the dataframe, and make it works for all the dictionaries
  columns=get_all_columns(alldicts,exclusions)
  #Set up dataframe with required columns
  df=pd.DataFrame(columns=columns)
  #Add the data to the dataframe
  for d in alldicts:
    fd=flatdict(d,columns,infofile,exclusions)
    df=df.append(fd,ignore_index=True)
  return df
  
def do_collection(inputfiles,outfpath,exclusions):
  """Generate a pandas dataframe from the specified input files, written to the output file
  Arguments:
    inputfiles = iterable of input files (yaml format) to read.
      This should be the complete paths, as strings.
    outfpath = path to the output file (pickle format), as a string.
    exclusions = list of keys in dictionary to ignore
  No return value.
  Output file is generated.
  """
  #Make sure the parent folder exists
  os.makedirs(osp.split(outfpath)[0],exist_ok=True)
  alldicts = [useful.readyaml(fp) for fp in inputfiles]
  df = dicts_to_dataframe(alldicts,exclusions)
  df.to_pickle(outfpath)
  return
