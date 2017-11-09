#Doit file for model runs

#Discussion
#This file is intended to automate the process from mesh generation, through running the solver, through post-processing.
#As such, it relies on scripts specific to each of those to generate the necessary tasks.
#The runs requested in control.yaml are used.

#Standard library
import os
import os.path as osp

#Site packages

#Local
from folderstructure import *
import useful
import tasks_mesh
import tasks_solver
import tasks_postproc

#Constants
controlfile=osp.join(datafolder,'control.yaml')

#Read in all the models and meshes
model_infiles=useful.readyaml(controlfile)
modelparams_filelist=[osp.join(params_model_folder,fn) for fn in model_infiles]
allmodels,modelfiles,allmeshes,meshfiles=tasks_solver.GetAllModelsAndMeshes(modelparams_filelist)

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
# def task_collect():
#   for fname in model_infiles:
#     basename = osp.splitext(fname)[0]
#     yield tasks_postproc.collection(basename,allmodels.values())

#Post-processing tasks
def task_postproc():
  return tasks_postproc.postproc_task_generator(model_infiles,allmodels)
