import os.path as osp

srcfolder=osp.split(__file__)[0]

#mesh
meshfolder=osp.join(srcfolder,'mesh')
geofolder=osp.join(meshfolder,'geo')
mshfolder=osp.join(meshfolder,'msh')
xmlfolder=osp.join(meshfolder,'xml')
outfolder=osp.join(meshfolder,'gmsh_out')

#solver
solverfolder=osp.join(srcfolder,'solver')

#solutions
solndir=osp.join(srcfolder,'solutions')