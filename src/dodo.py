#Doit file for model runs

#Discussion
#This file is intended to automate the process from mesh generation, through running the solver, through post-processing.
#As such, it relies on scripts specific to each of those to generate the necessary tasks.
#At present, the runs requested in control.yaml are used.
#This can be a symlink to a yaml file in the params folder.

#Standard library
from argparse import Namespace
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

#Read in the sequence of analysis runs
runs=AnalysisRunSet.all_from_yaml('control.yaml')

def consolidate(runs,yamlname,entryname,entrytype,otherprops=None):
  """Get list of all things of a certain type from all the AnalysisRunSet objects
  Explanation:
    control.yaml contains multiple yaml documents, each of which contains parameters with other yaml documents as their values
    You specify which yaml document you want to consolidate with argument 'yamlname',
    which must be the name of an AnalysisRunSet attribute.
    The yaml files thus specified also contain multiple documents.
    Each such document is here called an "entry".
    Argument 'entryname' specifies which of the entry's attributes specifies its name.
    Argument entrytype specifies the type of thing each entry actually is.
  Inputs:
    runs = iterator over AnalysisRunSet objects
    yamlname = attribute name in AnalysisRunSet providing the base name for the relevant yaml files
    entryname = attribute name in the entries themselves providing their names (because each yaml file can have more than one)
    entrytype = class to be used for each entry
    otherprops = name of any other properties from runs to be copied directly into each entry, as a sequence of attribute names (optional)
  Returns:
    entrylist = list of parameter objects of the given type"""
  entries_byname={}
  yaml_from_name={}
  for rd in runs:
    yamlfile=rd[yamlname]
    entry_iter=entrytype.all_from_yaml(yamlfile)
    for entry in entry_iter:
      if entry is not None:
        if otherprops is not None:
          for k in otherprops:
            setattr(entry,k,getattr(rd,k)
        obj=Namespace(**entry)
        objname=entry[entryname]
        assert objname not in entries_byname, "Duplicate %s name: %s in both %s and %s"%(entrytype.__name__,objname,yaml_from_name[objname],yamlfile)
        yaml_from_name[objname]=yamlfile
        entries_byname[objname]=obj
  entrylist=[x for x in entries_byname.values()]
  return entrylist

#Mesh tasks
def task_make_mesh():
  #Get list of all meshes to generate
  meshruns=consolidate(runs,'meshparams','meshname',buildgeom.MeshParameters)
  #Set up tasks for each mesh
  for params in meshruns:
    yield tasks_mesh.create_geo(params)
    yield tasks_mesh.create_msh(params)
    yield tasks_mesh.create_xml(params)

#Solver tasks
def task_solve():
  #Get list of all models to solve
  modelruns=consolidate(runs,'modelparams','modelname',solver_general.SolverParams,['meshparams'])
  #Set up tasks for each model
  for params in modelruns:
    yield tasks_solver.dosolve(params)
    
#Result collection tasks
def task_collect():
  #Get list of all results to collect
  for rd in runs:
    yield tasks_postproc.collection(rd)

