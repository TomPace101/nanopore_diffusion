
#Standard library
import os
import os.path as osp

#Site packages
from doit.tools import config_changed

#Local
import useful
import buildgeom

#Folders
topfolder=osp.abspath('mesh')
geofolder=osp.join(topfolder,'geo')
mshfolder=osp.join(topfolder,'msh')
xmlfolder=osp.join(topfolder,'xml')
outfolder=osp.join(topfolder,'gmsh_out')

def create_geo(params):
  paramdef=params.__dict__ #Ugly, but sometimes we need the dictionary version instead of an object
  geomyaml=osp.join(topfolder,params.lattice+'.yaml') #geometry definition yaml file
  geomdef=useful.readyaml(geomyaml) #geometry definition dictionary
  geofile=osp.join(geofolder,params.meshname+'.geo') #.geo file
  filedeps=[osp.join(topfolder,x) for x in ['buildgeom.py',geomyaml,geomdef['tmplfile'],'common.geo.jinja2'] ]
  tdef = {'name':geofile,
          'file_dep':filedeps,
          'uptodate':[config_changed(paramdef)],
          'targets':[geofile],
          'actions':[(buildgeom.write_one_geo,(geomdef,paramdef,geofile))]}
  return tdef

def create_msh(params):
  stem=params.meshname
  geofile=osp.join(geofolder,'%s.geo'%stem)
  mshfile=osp.join(mshfolder,'%s.msh'%stem)
  outfile=osp.join(outfolder,'%s.txt'%stem)
  tdef = {'name':mshfile,
          'file_dep':[geofile],
          'targets':[mshfile],
          'actions':['gmsh -0 -o %s %s >%s'%(mshfile,geofile,outfile)]}
  return tdef

def create_xml(params):
  stem=params.meshname
  mshfile=osp.join(mshfolder,'%s.msh'%stem)
  xmlfile=osp.join(xmlfolder,'%s.xml'%stem)
  tdef = {'name':xmlfile,
          'file_dep':[mshfile],
          'targets':[xmlfile],
          'actions':['dolfin-convert %s %s'%(mshfile,xmlfile)]}
  return tdef

