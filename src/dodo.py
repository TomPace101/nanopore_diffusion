#Doit file for model runs

#Discussion
#This file is intended to automate the process from mesh generation, through running the solver, through post-processing.
#As such, it relies on scripts specific to each of those to generate the necessary tasks.
#The runs requested in control.yaml are used.

#Standard library
import os
import os.path as osp
import sys

#Site packages

#Local
from folderstructure import *
import useful
import tasks_mesh
import tasks_solver
import tasks_postproc
import buildgeom
import solver_general

#Constants
controlfile='control.yaml'

def consolidate(entry_filelist,infolder,entrytype,nameattribute):
  """Load all objects from a list of files
  For each file in the list,
  this will load every document the file contains as an object of the specified type,
  and store all the objects from all the files in a single dictionary.
  Inputs:
    entry_filelist = list (or other iterable) of file names
    infolder = path to folder containing the files whose names are in entry_filelist
    entrytype = class to be used for each object
    nameattribute = attribute name in the entries providing their names
      These names are the keys in the output dictionaries.
  Returns:
    entries_byname = dictionary mapping an entry's name to the entry itself
    files_byname = dictionary mapping an entry's name to the file it came from."""
  entries_byname={}
  files_byname={}
  #For each listed file
  for filename in entry_filelist:
    #Get an iterable over the objects defined in this file
    infpath=osp.join(infolder,filename)
    entry_iter=entrytype.all_from_yaml(infpath)
    #For each object in this file
    for entry in entry_iter:
      if entry is not None:
        #Get its name
        objname=getattr(entry,nameattribute)
        #Do we already have one by that name?
        assert objname not in entries_byname, "Duplicate %s name: %s in both %s and %s"%(entrytype.__name__,objname,files_byname[objname],filename)
        #Add entry to dictionaries
        files_byname[objname]=infpath
        entries_byname[objname]=entry
  return entries_byname, files_byname

def GetAllModelsAndMeshes(controlfile):
  """Read in all the models and meshes
  Arguments:
    controlfile = path to file containing list of all model parameters files to read
  Return values:
    allmodels = Dictionary of all ModelParameters objects, by modelname
    modelfiles = Dictionary of yaml file for ModelParameters objects, by modelname
    allmeshses = Dictionary of all MeshParameters objects, by meshname
    meshfiles = Dictionary of yaml files for MeshParameters objects, by meshname"""
  #Read in the list of all model parameters files
  modelparams_filelist = AnalysisRunSet.readyaml(controlfile)
  #Get all the models from all the model parameter files
  allmodels, modelfiles = consolidate(modelparams_filelist,params_model_folder,solver_general.ModelParameters,'modelname')
  #Get a list of all mesh parameters files from the models
  meshparams_filelist = []
  for modelparams in allmodels.values():
    if not modelparams.meshparamsfile in meshparams_filelist:
      meshparams_filelist.append(modelparams.meshparams_filelist)
  #Get all the meshes from all the mesh parameter files
  allmeshes, meshfiles = consolidate(meshparams_filelist,params_mesh_folder,buildgeom.MeshParameters,'meshname')
  return allmodels,modelfiles,allmeshes,meshfiles

#Read in all the models and meshes
allmodels,modelfiles,allmeshes,meshfiles=GetAllModelsAndMeshes(controlfile)

#Mesh tasks
def task_make_mesh():
  #Set up tasks for each mesh
  for meshparams in allmeshes.values():
    yield tasks_mesh.create_geo(meshparams)
    yield tasks_mesh.create_msh(meshparams)
    yield tasks_mesh.create_xml(meshparams)

#Solver tasks
def task_solve():
  #Set up tasks for each model
  for modelparams in allmodels.values():
    meshparams=allmeshes[modelparams.meshname]
    yield tasks_solver.dosolve(modelparams,meshparams)
    
#Result collection tasks
def ignored_task_collect(): ##TODO: re-enable once working
  #Get list of all results to collect
  for r in runs:
    yield tasks_postproc.collection(r)

