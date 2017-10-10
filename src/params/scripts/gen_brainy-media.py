#Generate mesh/fig_Deff_fvf.yaml

#TODO: there's probably a lot that could be abstracted here.
#TODO: should this be included in doit somehow? (the whole parameter generation step, that is)

import itertools
import os.path as osp
import sys

import yaml

sys.path.append(osp.abspath('.'))
from folderstructure import *

hrhash='brainy-media'
cellsizes=[(12, 12), (24,6), (36,4), (5,5), (10,5)]
tmlist = [75, 25]
Hvals = {75:50, 25:25}
radstep = 0.5

dirichlet_pairs=[(5.0,1.0),(10.0,1.0)]

#The mesh parameters yaml file
meshfile=osp.join(params_mesh_folder,hrhash+'.yaml')
meshdocs=[]
idnum=1
for LxLy, tm in itertools.product(cellsizes,tmlist):
  assert idnum < 1000, "Insufficient number of digits provided"
  Rnow=radstep
  Lx, Ly = LxLy
  while Rnow < min(Lx, Ly):
    doc={'lattice':'body-cen2',
        'meshname':'%s_mesh_%03d'%(hrhash,idnum),
        'mscale': 1.0,
        'mcarh': 1.0,
        'mcarl': 5.0,
        'Lx': Lx,
        'Ly': Ly,
        'R': Rnow,
        'H': Hvals[tm],
        'tm': tm}
    meshdocs.append(doc)
    idnum += 1
    Rnow += radstep
with open(meshfile,'w') as fp:
  yaml.dump_all(meshdocs,fp)

meshname_list=[m['meshname'] for m in meshdocs]

#The model parameters file
modelfile=osp.join(params_model_folder,hrhash+'.yaml')
modeldocs=[]
idnum=1
for meshname,bcvals in itertools.product(meshname_list,dirichlet_pairs):
  assert idnum < 1000, "Insufficient number of digits provided"
  topval,baseval= bcvals
  doc={'modelname':'%s_model_%03d'%(hrhash,idnum),
        'meshname':meshname,
        'equation': 'fickian_unhomog',
        'topsurf': 4,
        'basesurf': 1,
        'topval': topval,
        'baseval': baseval,
        'fluxsurf': 12,
        'fluxsign': '-',
        'sample_spacing': 0.1}
  modeldocs.append(doc)
  idnum += 1
with open(modelfile,'w') as fp:
  yaml.dump_all(modeldocs,fp)

print('Generated %d meshes and %d models.'%(len(meshdocs),len(modeldocs)))