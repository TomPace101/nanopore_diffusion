#Functions used my various solvers

#Standard library
import os.path as osp
import sys

#Site packages

#Local
sys.path.append(osp.abspath('..'))
import useful
from folderstructure import *

def List_Mesh_Input_Files(params):
  mesh_xml=osp.join(xmlfolder,params.meshname+'.xml')
  surface_xml=osp.join(xmlfolder,params.meshname+'_facet_region.xml')
  volume_xml=osp.join(xmlfolder,params.meshname+'_physical_region.xml')
  return mesh_xml, surface_xml, volume_xml

