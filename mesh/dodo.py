
#TODO
#-the tasks are all defined before any are executed, so a new .geo file created won't be noticed.
#   you end up having to run doit 3 times to get everything updated
#-See "things that should probably come from somewhere else" below
#-related to these issues, we need a way to generate and name parametric variations automatically.
#   that means not reading the parameters from a single file like we do know
#   at the same time, I do want to have that option available.

#Standard library
import os
import os.path as osp

#Site packages
from doit.tools import config_changed

#Local
import useful
import buildgeom

#Folders
geofolder='geo'
mshfolder='msh'
xmlfolder='xml'
outfolder='gmsh_out'

#Things that should probably come from somewhere else
stemname = "test"
variants=['body','face']
paramfile='params.yaml'

def get_stemlist(folder,ext):
  "return the stem names of files with the given extension in the given folder (extension should be lowercase and begin with dot)"
  allfiles=os.listdir(folder)
  splitfiles=[osp.splitext(f) for f in allfiles]
  stemlist=[f for f,x in splitfiles if x.lower()==ext]
  return stemlist

def task_make_geo():
  paramdef=useful.readyaml(paramfile)
  for v in variants:
    geomyaml='%s-centered.yaml'%v
    geomdef=useful.readyaml(geomyaml)
    geofile=osp.join(geofolder,'%s-%s.geo'%(stemname,v))
    configdict={**paramdef, **geomdef} #python 3.5 and above only; for older versions do copy then update
    tdef={'name':'%s-%s'%(stemname,v),
          'file_dep':['buildgeom.py'],
          'uptodate':[config_changed(configdict)],
          'targets':[geofile],
          'actions':[(buildgeom.write_one_geo,(geomdef,paramdef,geofile))]}
    yield tdef

def task_make_msh():
  stemlist=get_stemlist(geofolder,'.geo')
  for stem in stemlist:
    geofile=osp.join(geofolder,'%s.geo'%stem)
    mshfile=osp.join(mshfolder,'%s.msh'%stem)
    outfile=osp.join(outfolder,'%s.txt'%stem)
    tdef = {'name':stem,
            'file_dep':[geofile],
            'targets':[mshfile],
            'actions':['gmsh -0 -o %s %s >%s'%(mshfile,geofile,outfile)]}
    yield tdef

def task_make_xml():
  stemlist=get_stemlist(mshfolder,'.msh')
  for stem in stemlist:
    mshfile=osp.join(mshfolder,'%s.msh'%stem)
    xmlfile=osp.join(xmlfolder,'%s.xml'%stem)
    tdef = {'name':stem,
            'file_dep':[mshfile],
            'targets':[xmlfile],
            'actions':['dolfin-convert %s %s'%(mshfile,xmlfile)]}
    yield tdef

