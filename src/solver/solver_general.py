#Functions used my various solvers

#Standard library
import os.path as osp
import sys

#Site packages

#Local
from folderstructure import *
import useful

class SolverParams(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store generic solver parameters
  Attributes:
    modelname = stem name for output files
    meshname = stem name for mesh files
    equation = name of equation to be solved
    boundaryconditions = parameters specifying boundary conditions
      The parameters specified are specific to the euqation being solved
    dataextraction = parameters for data to extract from the solution"""
  __slots__=('modelname','meshname','equation','boundaryconditions','dataextraction')

def List_Mesh_Input_Files(params):
  mesh_xml=osp.join(xmlfolder,params.meshname+'.xml')
  surface_xml=osp.join(xmlfolder,params.meshname+'_facet_region.xml')
  volume_xml=osp.join(xmlfolder,params.meshname+'_physical_region.xml')
  return mesh_xml, surface_xml, volume_xml

