#Doit file for model runs

#TODO
#This still seems too complicated.

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
import useful
from mesh import tasks_mesh

#Read the yaml document
runs=useful.readyaml_multidoc('control.yaml')

#Get list of all meshes to generate
meshes_byname={}
yaml_from_meshname={}
for rd in runs:
  meshyaml=rd['meshparams']
  meshdict_list=useful.readyaml_multidoc(meshyaml)
  for meshdict in meshdict_list:
    if meshdict is not None:
      meshobj=Namespace(**meshdict)
      meshname=meshdict['meshname']
      assert meshname not in meshes_byname, "Duplicate mesh name: %s in both %s and %s"%(meshname,yaml_from_meshname[meshname],meshyaml)
      yaml_from_meshname[meshname]=meshyaml
      meshes_byname[meshname]=meshobj
meshruns=[m for m in meshes_byname.values()]

#Mesh tasks
def task_make_mesh():
  for params in meshruns:
    yield tasks_mesh.create_geo(params)
    yield tasks_mesh.create_msh(params)
    yield tasks_mesh.create_xml(params)

