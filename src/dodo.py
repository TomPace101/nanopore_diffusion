
#TODO

#Discussion
#This file is intended to automate the process from mesh generation, through running the solver, through post-processing.
#As such, it relies on scripts specific to each of those to generate the necessary tasks.

#Standard library
import os
import os.path as osp

#Site packages

#Local
sys.path.append(osp.abspath('.'))
import useful
sys.path.append(osp.abspath('./mesh'))
import tasks_mesh

#Read the yaml document
runs=readyaml_multidoc('control.yaml')

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
