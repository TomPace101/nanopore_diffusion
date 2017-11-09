
#Standard library
import os.path as osp
import sys

srcfolder=osp.abspath(osp.split(__file__)[0])
datafolder=osp.join(osp.split(srcfolder)[0],'data')

#params
paramsfolder=osp.join(datafolder,'params')
params_mesh_folder=osp.join(paramsfolder,'mesh')
params_model_folder=osp.join(paramsfolder,'model')
params_postproc_folder=osp.join(paramsfolder,'postproc')

#mesh
meshfolder=osp.join(datafolder,'mesh')
geomdef_folder=osp.join(meshfolder,'geomdef')
geotemplates_folder=osp.join(meshfolder,'templates')
geofolder=osp.join(meshfolder,'geo')
mshfolder=osp.join(meshfolder,'msh')
xmlfolder=osp.join(meshfolder,'xml')
gmsh_outfolder=osp.join(meshfolder,'gmsh_out')

#solutions
solnfolder=osp.join(datafolder,'solutions')

#post-processing
postprocfolder=osp.join(datafolder,'postproc')

#add python code folder(s) to path
if not srcfolder in sys.path:
  sys.path.append(srcfolder)

