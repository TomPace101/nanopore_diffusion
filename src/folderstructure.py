
#Standard library
import os
import os.path as osp
import sys

#Locate source folders
if 'SRCLOC' in os.environ.keys():
  srcfolder=osp.normpath(osp.abspath(os.environ['SRCLOC']))
else:
  srcfolder=osp.abspath(osp.split(__file__)[0])

#Locate data folder
if 'DATALOC' in os.environ.keys():
  datafolder=osp.normpath(osp.abspath(os.environ['DATALOC']))
else:
  datafolder=osp.join(osp.split(srcfolder)[0],'data')

#subfolders of src
custom_modules_folder=osp.join(srcfolder,'customizations')
solver_modules_folder=osp.join(srcfolder,'solvers')

#params
paramsfolder=osp.join(datafolder,'params')
params_mesh_folder=osp.join(paramsfolder,'mesh')
params_model_folder=osp.join(paramsfolder,'model')
params_postproc_folder=osp.join(paramsfolder,'postproc')
params_paramgen_folder=osp.join(paramsfolder,'paramgen')

#mesh
meshfolder=osp.join(datafolder,'mesh')
geomdef_folder=osp.join(meshfolder,'geomdef')
geotemplates_folder=osp.join(meshfolder,'templates')
geofolder=osp.join(meshfolder,'geo')
mshfolder=osp.join(meshfolder,'msh')
xmlfolder=osp.join(meshfolder,'xml')
gmsh_outfolder=osp.join(meshfolder,'gmsh_out')
paramlocs_outfolder=osp.join(meshfolder,'paramlocs')

#solutions
solnfolder=osp.join(datafolder,'solutions')

#post-processing
postprocfolder=osp.join(datafolder,'postproc')

#parameter generation
pgtemplates_folder=osp.join(datafolder,'paramgen_tmpl')

#add python code folder(s) to path
if not srcfolder in sys.path:
  sys.path.append(srcfolder)
if not custom_modules_folder in sys.path:
  sys.path.append(custom_modules_folder)
if not solver_modules_folder in sys.path:
  sys.path.append(solver_modules_folder)
