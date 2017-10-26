
#Standard library
import os
import os.path as osp

#Site packages
from doit.tools import config_changed

#Local
from folderstructure import *
import useful
import buildgeom

def assure_dir(fpath):
  "Create folder for a given path, if it does not exist"
  if not osp.isdir(osp.dirname(fpath)):
    os.makedirs(osp.dirname(fpath))
  return

def assure_multi_dirs(*args):
  "Call assure_dir for each path in a list of arguments"
  for fpath in args:
    assure_dir(fpath)
  return

def create_geo(params):
  geomyaml, geomdef, geofile = buildgeom.process_mesh_params(params)
  filedeps=[osp.join(meshfolder,'buildgeom.py'),
            osp.join(geomdef_folder,geomyaml),
            osp.join(geotemplates_folder,geomdef.tmplfile),
            osp.join(geotemplates_folder,'common.geo.jinja2')]
  tdef = {'name':params.meshname+'.geo',
          'file_dep':filedeps,
          'uptodate':[config_changed(params.to_dict())],
          'targets':[geofile],
          'actions':[(buildgeom.write_one_geo,(geomdef,params,geofile))]}
  return tdef

def create_msh(params):
  stem=params.meshname
  geofile=osp.join(geofolder,params.basename,'%s.geo'%stem)
  mshfile=osp.join(mshfolder,params.basename,'%s.msh'%stem)
  outfile=osp.join(gmsh_outfolder,params.basename,'%s.txt'%stem)
  tdef = {'name':stem+'.msh',
          'file_dep':[geofile],
          'targets':[mshfile],
          'actions':[(assure_multi_dirs,(mshfile,outfile)),'gmsh -0 -o %s %s >%s'%(mshfile,geofile,outfile)]}
  return tdef

def create_xml(params):
  stem=params.meshname
  mshfile=osp.join(mshfolder,params.basename,'%s.msh'%stem)
  xmlfile=osp.join(xmlfolder,params.basename,'%s.xml'%stem)
  tdef = {'name':stem+'.xml',
          'file_dep':[mshfile],
          'targets':[xmlfile],
          'actions':[(assure_dir,(xmlfile,)), 'dolfin-convert %s %s'%(mshfile,xmlfile)]}
  return tdef

