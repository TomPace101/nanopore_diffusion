
#Standard library
import os
import os.path as osp
import sys

#Site packages
from doit.tools import config_changed

#Local
sys.path.append('..')
from folderstructure import *
import useful
import buildgeom

def create_geo(params):
  geomyaml, geomdef, geofile = buildgeom.process_mesh_params(params)
  filedeps=[osp.join(meshfolder,'buildgeom.py'),
            osp.join(geomdef_folder,geomyaml),
            osp.join(geotemplates_folder,geomdef.tmplfile),
            osp.join(geotemplates_folder,'common.geo.jinja2')]
  tdef = {'name':geofile,
          'file_dep':filedeps,
          'uptodate':[config_changed(params.to_dict())],
          'targets':[geofile],
          'actions':[(buildgeom.write_one_geo,(geomdef,params,geofile))]}
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

