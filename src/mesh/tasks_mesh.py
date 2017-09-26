
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

def create_geo(paramdef):
  params=useful.dict_to_nt(paramdef) #For convenience
  geomyaml=osp.join(topfolder,params.lattice+'.yaml') #geometry definition yaml file
  geomdef=useful.readyaml(geomyaml) #geometry definition dictionary
  geofile=osp.join(geofolder,params.basename+'.geo') #.geo file
  configdict={**paramdef, **geomdef} #python 3.5 and above only; for older versions do copy then update
  tdef = {'name':params.basename,
          'file_dep':['buildgeom.py',geomdef['tmplfile'],osp.join(topfolder,'common.geo.jinja2')],
          'uptodate':[config_changed(configdict)],
          'targets':[geofile],
          'actions':[(buildgeom.write_one_geo,(geomdef,paramdef,geofile))]}
  return tdef

def create_msh(paramdef):
  stem=paramdef['basename']
  geofile=osp.join(geofolder,'%s.geo'%stem)
  mshfile=osp.join(mshfolder,'%s.msh'%stem)
  outfile=osp.join(outfolder,'%s.txt'%stem)
  tdef = {'name':stem,
          'file_dep':[geofile],
          'targets':[mshfile],
          'actions':['gmsh -0 -o %s %s >%s'%(mshfile,geofile,outfile)]}
  return tdef

def create_xml(paramdef):
  stem=paramdef['basename']
  mshfile=osp.join(mshfolder,'%s.msh'%)
  xmlfile=osp.join(xmlfolder,'%s.xml'%stem)
  tdef = {'name':stem,
          'file_dep':[mshfile],
          'targets':[xmlfile],
          'actions':['dolfin-convert %s %s'%(mshfile,xmlfile)]}
  return tdef

