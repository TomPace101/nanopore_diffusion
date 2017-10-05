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
from solver import tasks_solver

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
      meshname=meshobj.meshname
      assert meshname not in meshes_byname, "Duplicate mesh name: %s in both %s and %s"%(meshname,yaml_from_meshname[meshname],meshyaml)
      yaml_from_meshname[meshname]=meshyaml
      meshes_byname[meshname]=meshobj
meshruns=[m for m in meshes_byname.values()]

#Get list of all models to solve
bcs_byname={}
yaml_from_bcname={}
for rd in runs:
  bcyaml=rd['bcparams']
  bcdict_list=useful.readyaml_multidoc(bcyaml)
  for bcdict in bcdict_list:
    if bcdict is not None:
      bcobj=Namespace(**bcdict)
      bcname=bcobj.bcname
      assert bcname not in bcs_byname, "Duplicate boundary condition name: %s in both %s and %s"%(bcname,yaml_from_bcname[bcname],bcyaml)
      yaml_from_bcname[bcname]=bcyaml
      bcs_byname[bcname]=bcobj
bcruns=[bc for bc in bcs_byname.values()]

#Mesh tasks
def task_make_mesh():
  for params in meshruns:
    yield tasks_mesh.create_geo(params)
    yield tasks_mesh.create_msh(params)
    yield tasks_mesh.create_xml(params)

#Solver tasks
def task_solve():
  for params in bcruns:
    yield tasks_solver.dosolve(params)
