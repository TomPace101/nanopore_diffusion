#Functions used my various solvers

#Standard library
import os.path as osp
import sys

#Site packages

#Local
from folderstructure import *
import useful
import extraction_functions

class ModelParameters(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store generic solver parameters
  Attributes:
    modelname = stem name for output files
    meshname = stem name for mesh files
    equation = name of equation to be solved
    boundaryconditions = parameters specifying boundary conditions
      The parameters specified are specific to the euqation being solved
    dataextraction = parameters for data to extract from the solution"""
  __slots__=('modelname','meshname','equation','boundaryconditions','dataextraction')

def List_Mesh_Input_Files(meshname):
  mesh_xml=osp.join(xmlfolder,meshname+'.xml')
  surface_xml=osp.join(xmlfolder,meshname+'_facet_region.xml')
  volume_xml=osp.join(xmlfolder,meshname+'_physical_region.xml')
  return mesh_xml, surface_xml, volume_xml

def Create_Output(modelparams,meshparams,soln,cmdlist):
  """Process a sequence of data extraction commands on the given solution.
  Arguments:
    modelparams = ModelParameters object
    meshparams = buildgeom.MeshParameters object
    soln = the FEniCS solution to extract data from
    cmdlist = sequence of extraction commands
      Each command is a pair of extraction function names and keyword arguments.
  No return value.
  Output files are generated."""
  #Output location(s)
  outdir=osp.join(solnfolder,modelparams.modelname)
  if not osp.isdir(outdir):
    os.mkdir(outdir)

  #Initialize results dictionary
  results=modelparams.to_dict()
  results.update(meshparams.to_dict())

  #Process each command
  for cmd in cmdlist:
    #Function name and arguments
    funcname, kwargs = cmd
    #Call it
    extraction_functions.exfuncs[funcname](soln,results,outdir,**kwargs)
    
  #Write out the results file
  useful.writeyaml(results,osp.join(outdir,'results.yaml'))

  