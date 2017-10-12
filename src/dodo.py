#Doit file for model runs

#Discussion
#This file is intended to automate the process from mesh generation, through running the solver, through post-processing.
#As such, it relies on scripts specific to each of those to generate the necessary tasks.
#At present, the runs requested in control.yaml are used.
#This can be a symlink to a yaml file in the params folder.

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

class AnalysisRunSet(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store the data for completing a set of analysis runs
  Attributes:
    meshparams = yaml file containing buildgeom.MeshParameters documents
    modelparams = yaml file containing 
    master = """
  __slots__=('meshparams','modelparams','master')

def consolidate(runs,yamlfolder,yamlname,entrytype,nameattribute):
  """Get all things of a certain type from all the AnalysisRunSet objects
  Explanation:
    control.yaml contains multiple AnalysisRunSet objects, many of which specify other yaml files,
    which in turn contain multiple objects of various types.
    This function reads all of the yaml files for a certain type of object,
    and consolidates those objects into one dictionary.
  Inputs:
    runs = iterable over AnalysisRunSet objects
    yamlfolder = path to folder containing yaml files referenced by yamlname, as string
    yamlname = attribute name in AnalysisRunSet providing the base name for the relevant yaml files
    entrytype = class to be used for each object
    nameattribute = attribute name in the entries themselves providing their names
      These names are the keys in the output dictionary.
  Returns:
    entries_byname = dictionary mapping an entry's name to the entry itself
    yaml_from_name = dictionary mapping an entry's name to the yaml file it came from"""
  entries_byname={}
  yaml_from_name={}
  for rd in runs:
    yamlfile=getattr(rd,yamlname)+'.yaml'
    yamlfpath=osp.join(yamlfolder,yamlfile)
    entry_iter=entrytype.all_from_yaml(yamlfpath)
    for entry in entry_iter:
      if entry is not None:
        objname=getattr(entry,nameattribute)
        assert objname not in entries_byname, "Duplicate %s name: %s in both %s and %s"%(entrytype.__name__,objname,yaml_from_name[objname],yamlfile)
        yaml_from_name[objname]=yamlfpath
        entries_byname[objname]=entry
  return entries_byname, yaml_from_name

#Read in the sequence of analysis runs
runs=[r for r in AnalysisRunSet.all_from_yaml('control.yaml')]

#Get all the meshes and models
allmeshes,meshfiles=consolidate(runs,params_mesh_folder,'meshparams',buildgeom.MeshParameters,'meshname')
allmodels,modelfiles=consolidate(runs,params_model_folder,'modelparams',solver_general.ModelParameters,'modelname')

#Mesh tasks
def task_make_mesh():
  #Set up tasks for each mesh
  for params in allmeshes.values():
    yield tasks_mesh.create_geo(params)
    yield tasks_mesh.create_msh(params)
    yield tasks_mesh.create_xml(params)

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

