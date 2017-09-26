#Doit file for model runs

#TODO
#This still seems too complicated.

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
sys.path.append(osp.abspath('.'))
import useful
sys.path.append(osp.abspath('./mesh'))
import tasks_mesh

#Read the yaml document
runs=useful.readyaml_multidoc('control.yaml')

#Mesh tasks
def task_make_geo():
  for paramdef in runs:
    if paramdef is not None:
      yield tasks_mesh.create_geo(paramdef)

def task_make_msh():
  for paramdef in runs:
    if paramdef is not None:
      yield tasks_mesh.create_msh(paramdef)

def task_make_xml():
  for paramdef in runs:
    if paramdef is not None:
      yield tasks_mesh.create_xml(paramdef)
