#Standard library
import importlib
import os.path as osp
import sys

#Site packages

#Local
import folderstructure as FS
import common
import buildgeom
import solver_general

#Solver modules
solver_module_list=['solver_fickian_unhomog','solver_smol_unhomog']

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

def list_outputfiles(cmdlist):
  """Get a list of all the files generated by the data extraction commands.
  The info.yaml file is included as well.
  Arguments:
    cmdlist = list of data extraction commands,
      each command consists of pair (cmdname, arguments):
        cmdname = name of data extraction method of the solver class
        arguments = dictionary of all arguments needed by the extraction method
  Return:
    outfiles = list of generated output files (names only, not including their folder)"""
  #Currently, we assume all files can only come from the 'filename' argument
  filearg_list=['filename']
  outfiles=['info.yaml']
  for cmdname, arguments in cmdlist:
    #Check all possible arguments that could contain the name of an output file
    present_args=[n for n in filearg_list if n in arguments.keys()]
    outfiles.extend([arguments[n] for n in present_args])
  return outfiles

class ModelParameters(solver_general.ModelParametersBase):
  """Extend ModelParametersBase to allow for task generation and execution
  Attributes:
    solverclass = the class used to solve the model
    meshparams = instance of buildgeom.MeshParameters
    xmlfolder = path to xml input files for mesh
    mesh_xml, surface_xml, volume_xml = various mesh input xml files, full paths, as strings"""
  __slots__=('solverclass','meshparams','xmlfolder','mesh_xml','surface_xml','volume_xml','_more_inputfiles','_more_outputfiles')
  _folders={'meshparamsfile':FS.params_mesh_folder}
  #don't need sourcefile as input file due to config
  _inputfile_attrs=['mesh_xml', 'surface_xml', 'volume_xml']
  #Remember loaded meshes (parameters, not xml) so we aren't reading the same files over and over when solving multiple models
  loaded_meshfiles=[]
  loaded_meshes={}
  #Load the solver modules and map equation names to their solver classes
  #Mapping from ModelParameters.equation to the appropriate solver classes
  #Each module implementing a solver should define a solverclasses dictionary, so we just need to put them all together
  solverclasses={}
  for sm_name in solver_module_list:
    sm=importlib.import_module(sm_name)
    solverclasses.update(sm.solverclasses)
  
  def __init__(self,**kwd):
    #Initialization from base class
    super().__init__(**kwd)
    #Load the meshparameters file if not already loaded
    if not self.meshparamsfile in self.loaded_meshfiles:
      meshparams_gen=buildgeom.MeshParameters.all_from_yaml(self.full_path('meshparamsfile'))
      for mp in meshparams_gen:
        self.loaded_meshes[mp.meshname]=mp
      self.loaded_meshfiles.append(self.meshparamsfile)
    #Get the requested mesh
    self.meshparams=self.loaded_meshes[self.meshname]
    #Now get the location of the XML 
    #Get the solver class
    self.solverclass=self.solverclasses[self.equation]
    #Get code files
    self._more_inputfiles=[thisfile,common.__file__, sys.modules[self.solverclass.__module__].__file__]
    #Get XML files
    self.xmlfolder=osp.join(FS.xmlfolder,self.meshparams.basename)
    self.mesh_xml=osp.join(self.xmlfolder,self.meshparams.meshname+'.xml')
    self.surface_xml=osp.join(self.xmlfolder,self.meshparams.meshname+'_facet_region.xml')
    self.volume_xml=osp.join(self.xmlfolder,self.meshparams.meshname+'_physical_region.xml')
    #Get output files
    outfiles=list_outputfiles(getattr(self,'dataextraction',[]))
    self._more_outputfiles=[osp.join(self.outdir,f) for f in outfiles]

  def run(self):
    """Run the loaded model."""
    #Run the solver
    print(self.modelname)
    self.solverclass.complete(self,self.meshparams)


#Support command-line arguments
if __name__ == '__main__':
  program_description='Solve a diffusion equation with fenics'
  input_file_description='Path to file containing ModelParameters definitions'
  other_selection={'equation':ModelParameters.solverclasses.keys()}
  
  common.run_cmd_line(program_description,input_file_description,ModelParameters,other_selection=other_selection)
  #other_selection is needed so we only try to run models whose equation we have a solver for
