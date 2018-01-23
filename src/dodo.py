#Doit file for model runs

#Discussion
#This file is intended to automate the process from mesh generation, through running the solver, through post-processing.
#As such, it relies on scripts specific to each of those to generate the necessary tasks.
#The runs requested in control.yaml are used.

#Standard library
import os
import os.path as osp

#Site packages
from doit import create_after

#Local
import folderstructure as FS
import useful
import paramgen
import buildgeom
import geom_mk_msh
import geom_mk_xml
import solver_run

#Constants
controlfile=osp.join(FS.datafolder,'control.yaml')

#Read in the files for processing
infile_list=useful.readyaml(controlfile)


def generic_task_generator(folder,objtype,infile_list=infile_list):
  for infile in infile_list:
    infpath=osp.join(folder,infile)
    if osp.isfile(infpath):
      print("Loading %s from %s."%(objtype.__name__,infpath))
      allobj=objtype.all_from_yaml(infpath)
      for obj in allobj:
        yield obj.task_definition

#Parameter generation tasks
def task_paramgen():
  return generic_task_generator(FS.params_paramgen_folder,paramgen.ParameterGenerator)

#Mesh tasks
@create_after('paramgen')
def task_make_geo():
  return generic_task_generator(FS.params_mesh_folder,buildgeom.MeshParameters)

@create_after('paramgen')
def task_make_msh():
  return generic_task_generator(FS.params_mesh_folder,geom_mk_msh.GmshRunner)

@create_after('paramgen')
def task_make_xml():
  return generic_task_generator(FS.params_mesh_folder,geom_mk_xml.DolfinConvertRunner)

#Solver tasks
@create_after('paramgen')
def task_solve():
  return generic_task_generator(FS.params_model_folder,solver_run.ModelParameters)

##TODO
# def task_solve():
#   for fn in infile_list:
#     models=loadobjs(folderstructure.params_model_folder,fn,solver_general.ModelParameters)
#     for modelparams in models:
#       ##TODO
#       pass

#Post-processing tasks
##TODO




##------------------------------------------------------------------------------
##TODO: everything below here is old and should be deleted

##import tasks_mesh
##import tasks_solver
##import tasks_postproc


# def loadobjs(indir,fname,objtype):
#   infpath = osp.join(indir,fname)
#   if osp.isfile(infpath):
#     gen=objtype.all_from_yaml(infpath)
#   else:
#     gen=[]
#   return gen

# modelparams_filelist=[osp.join(params_model_folder,fn) for fn in infile_list]
# allmodels,modelfiles,allmeshes,meshfiles=tasks_solver.GetAllModelsAndMeshes(modelparams_filelist)
# 
# #Solver tasks
# def task_solve():
#   #Set up tasks for each model
#   for modelparams in allmodels.values():
#     meshparams=allmeshes[modelparams.meshname]
#     yield tasks_solver.dosolve(modelparams,meshparams)
#     
# #Result collection tasks
# # def task_collect():
# #   for fname in infile_list:
# #     basename = osp.splitext(fname)[0]
# #     yield tasks_postproc.collection(basename,allmodels.values())
# 
# #Post-processing tasks
# def task_postproc():
#   return tasks_postproc.postproc_task_generator(infile_list,allmodels)
