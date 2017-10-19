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
import solver_general

#Constants
controlfile='control.yaml'


#Read in all the models and meshes
modelparams_filelist=useful.readyaml(controlfile)
allmodels,modelfiles,allmeshes,meshfiles=solver_general.GetAllModelsAndMeshes(modelparams_filelist)

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

