
#Standard library
import os.path as osp
import sys

srcfolder=osp.split(__file__)[0]

#params
paramsfolder=osp.join(srcfolder,'params')
params_mesh_folder=osp.join(paramsfolder,'mesh')
params_model_folder=osp.join(paramsfolder,'model')
params_control_folder=osp.join(paramsfolder,'control')

#mesh
meshfolder=osp.join(srcfolder,'mesh')
geomdef_folder=osp.join(meshfolder,'geomdef')
geotemplates_folder=osp.join(meshfolder,'templates')
geofolder=osp.join(meshfolder,'geo')
mshfolder=osp.join(meshfolder,'msh')
xmlfolder=osp.join(meshfolder,'xml')
outfolder=osp.join(meshfolder,'gmsh_out')

#solver
solverfolder=osp.join(srcfolder,'solver')

#solutions
solnfolder=osp.join(srcfolder,'solutions')

#postproc
postprocfolder=osp.join(srcfolder,'postproc')

#add python code folders to path
sys.path.append(srcfolder)
sys.path.append(meshfolder)
sys.path.append(solverfolder)
sys.path.append(postprocfolder)

