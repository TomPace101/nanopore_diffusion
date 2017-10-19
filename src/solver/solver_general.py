#Functions used my various solvers

#Standard library
import os
import os.path as osp
import sys

#Site packages
from fenics import *
import numpy as np

#Local
from folderstructure import *
import useful
import plotdata
import buildgeom


class ModelParameters(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store generic solver parameters
  Attributes:
    modelname = stem name for output files
    meshparamsfile = name of yaml file containing the mesh defined by meshname (include .yaml extension)
    meshname = stem name for mesh files
    equation = name of equation to be solved
    properties = dictionary of property values
    boundaryconditions = parameters specifying boundary conditions
      The parameters specified are specific to the euqation being solved
    dataextraction = parameters for data to extract from the solution"""
  __slots__=('modelname','meshparamsfile','meshname','equation','properties','boundaryconditions','dataextraction')

def List_Mesh_Input_Files(meshname):
  mesh_xml=osp.join(xmlfolder,meshname+'.xml')
  surface_xml=osp.join(xmlfolder,meshname+'_facet_region.xml')
  volume_xml=osp.join(xmlfolder,meshname+'_physical_region.xml')
  return mesh_xml, surface_xml, volume_xml

def consolidate(entry_filelist,infolder,entrytype,nameattribute):
  """Load all objects from a list of files
  For each file in the list,
  this will load every document the file contains as an object of the specified type,
  and store all the objects from all the files in a single dictionary.
  Inputs:
    entry_filelist = list (or other iterable) of file names
    infolder = path to folder containing the files whose names are in entry_filelist
    entrytype = class to be used for each object
    nameattribute = attribute name in the entries providing their names
      These names are the keys in the output dictionaries.
  Returns:
    entries_byname = dictionary mapping an entry's name to the entry itself
    files_byname = dictionary mapping an entry's name to the file it came from."""
  entries_byname={}
  files_byname={}
  #For each listed file
  for filename in entry_filelist:
    #Get an iterable over the objects defined in this file
    infpath=osp.join(infolder,filename)
    entry_iter=entrytype.all_from_yaml(infpath)
    #For each object in this file
    for entry in entry_iter:
      if entry is not None:
        #Get its name
        objname=getattr(entry,nameattribute)
        #Do we already have one by that name?
        assert objname not in entries_byname, "Duplicate %s name: %s in both %s and %s"%(entrytype.__name__,objname,files_byname[objname],filename)
        #Add entry to dictionaries
        files_byname[objname]=infpath
        entries_byname[objname]=entry
  return entries_byname, files_byname

def ListMeshParamsFiles(modelparams_list):
  """Get a list of all mesh parameters files from the models
  Arguments:
    modelparams_list: list (or other iterable) of ModelParameters objects
  Returns:
    meshparams_filelist = list of mesh parameters filenames (not including folder path)"""
  meshparams_filelist = []
  for modelparams in modelparams_list:
    if not modelparams.meshparamsfile in meshparams_filelist:
      meshparams_filelist.append(modelparams.meshparamsfile)
  return meshparams_filelist

def GetAllModelsAndMeshes(modelparams_filelist):
  """Read in all the models and meshes
  Arguments:
    modelparams_filelist = list of model parameters files
  Return values:
    allmodels = Dictionary of all ModelParameters objects, by modelname
    modelfiles = Dictionary of yaml file for ModelParameters objects, by modelname
    allmeshses = Dictionary of all MeshParameters objects, by meshname
    meshfiles = Dictionary of yaml files for MeshParameters objects, by meshname"""
  #Get all the models from all the model parameter files
  allmodels, modelfiles = consolidate(modelparams_filelist,params_model_folder,ModelParameters,'modelname')
  #Get a list of all mesh parameters files from the models
  meshparams_filelist = ListMeshParamsFiles(allmodels.values())
  #Get all the meshes from all the mesh parameter files
  allmeshes, meshfiles = consolidate(meshparams_filelist,params_mesh_folder,buildgeom.MeshParameters,'meshname')
  return allmodels,modelfiles,allmeshes,meshfiles

class GenericSolver:
  """A generic solver, to be subclassed by solvers for the specific equations
  Subclasses should, at a minimum, implement a "solve" method to generate the solution
  and other data needed by their data extraction functions.
  Subclasses may choose to override the extraction functions provided here.
  Attributes:
    modelparams = ModelParameters instance
    meshparams = buildgeom.MeshParameters instance
    mesh = FEniCS Mesh
    surfaces = FEniCS MeshFunction of Physical Surface nubmer
    volumes = FEniCS MeshFunction of Physical Volume number
  """
  def __init__(self,modelparams,meshparams):
    """Initialize the solver by loading the Mesh and MeshFunctions.
    Arguments:
      modelparams = ModelParameters instance
      meshparams = buildgeom.MeshParameters instance"""
    #Store defining ParameterSet objects
    self.modelparams=modelparams
    self.meshparams=meshparams
      
    #Mesh input files
    mesh_xml, surface_xml, volume_xml = List_Mesh_Input_Files(modelparams.meshname)

    #Load mesh and meshfunctions
    self.mesh=Mesh(mesh_xml)
    self.surfaces=MeshFunction("size_t", mesh, surface_xml) #Mesh Function of Physical Surface number
    self.volumes=MeshFunction("size_t", mesh, volume_xml) #Mesh function of Physical Volume number

  def complete(self):
    "Convenience function to solve the model and generate all the requested output."
    self.solve()
    self.create_output()
    return

  def solve(self):
    "Method to be overridden by derived classes"
    raise NotImplementedError("%s did not override 'solve' method."%str(type(self)))

  def create_output(self):
    """Process the data extraction commands
    Adds the following attributes:
      outdir = path to output directory, as string
      results = dictionary of input and output values
    Output files are generated."""
    #Output location(s)
    self.outdir=osp.join(solnfolder,self.modelparams.modelname)
    if not osp.isdir(self.outdir):
      os.mkdir(self.outdir)

    #Initialize results dictionary
    self.results=self.modelparams.to_dict()
    self.results.update(self.meshparams.to_dict())

    #Process each command
    for cmd in self.modelparams.dataextraction:
      #Function name and arguments
      funcname, kwargs = cmd
      #Call it
      getattr(self,funcname)(**kwargs)
      extraction_functions.exfuncs[funcname](soln,results,outdir,**kwargs)
      
    #Write out the results file
    useful.writeyaml(self.results,osp.join(self.outdir,'results.yaml'))

    #Done
    return

  def solutionfield(filename):
    """Write solution field to VTK file
    Arguments:
      filename = name of output file, as string
        File will be created in the output directory (self.outdir)
    No attributes from other data extractions required.
    No new attributes.
    No return value.
    Output file is written."""
    vtk_file=File(osp.join(self.outdir,filename))
    vtk_file << self.soln
    return

  def fluxfield(filename): ##TODO: need D_bulk to get correct units
    """Flux as vector field (new attribute, and VTK file)
    Arguments:
      filename = name of output file, as string
        File will be created in the output directory (self.outdir)
    No attributes from other data extractions required.
    New attributes:
      flux = flux field as a MeshFunction
    No return value.
    Output file is written."""
    D_bulk=self.modelparams.properties['D_bulk']
    self.flux=project(Constant(-D_bulk)*grad(self.soln),self.V_vec)
    vtk_file=File(osp.join(self.outdir,filename))
    vtk_file << self.fluxfield
    return

  def fluxintegral(fluxsurf,fluxsign,name): ##TODO: store also quadrupled value for unit cell
    """Flux integral over specified surface
    Arguments:
      fluxsurf = physical surface number for flux measurement
      fluxsign = '+' or '-' to specify which diretion normal to the surface for flux calculation
      name = name for storage in the results dictionary
    Required attributes:
      flux = as calculated by fluxfield
    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    n=FacetNormal(self.mesh)
    dsi=Measure('interior_facet',domain=self.mesh,subdomain_data=self.surfaces)
    totflux=assemble(dot(self.flux,n(fluxsign))*dsi(fluxsurf))
    self.results[name]=totflux
    return

  def effective_diffusion(name,totflux_name):
    """Calculate effective diffusion constant
    Arguments:
      name = name for storage in the results dictionary
      totflux_name = name of previously calculated total flux in results dictionary
        This requires a prevous call to fluxintegral.
    No attributes from other data extractions required.
    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    quarter_area = self.meshparams.Lx * self.meshparams.Ly
    samples=[soln(0,0,zv) for zv in [self.meshparams.H, self.meshparams.H + self.meshparams.tm]]
    delta=samples[1]-samples[0]
    Deff=float(self.results[totflux_name]/quarter_area*self.meshparams.tm/delta)
    results[name]=Deff
    return

  def volfrac(soln,results,outdir,name):
    """Calculate free volume fraction
    Arguments:
      name = name for storage in the results dictionary
    No attributes from other data extractions required.
    No new attributes.
    New item added to results dictionary.
    No return value.
    No output files."""
    results[name]=np.pi*self.meshparams.R**2/(4*self.meshparams.Lx*self.meshparams.Ly)
    return

  def profile_centerline(spacing,filename,label):
    """Data for plot of concentration profile along centerline
    Arguments:
      spacing = distance between sampled points for line plots
      filename = name of output file, as string
        File will be created in the output directory (self.outdir)
      label = series label to assign, as string
    No attributes from other data extractions required.
    No new attributes.
    New item added to results dictionary.
    No return value.
    Output file is written."""
    ##TODO: replace c with soln
    zr=np.arange(0, 2*self.meshparams.H + self.meshparams.tm, spacing)
    zlist=[]
    vlist=[]
    for z in zr:
      zlist.append(z)
      tup=(0,0,z)
      vlist.append(soln(*tup))
    zarr=np.array(zlist)
    varr=np.array(vlist)
    meta=dict([(k,getattr(self.meshparams,k)) for k in ['Lx','Ly','R','tm','H']])
    meta.update(self.modelparams.boundaryconditions)
    pd=plotdata.PlotSeries(xvals=zarr,yvals=varr,label=label,metadata=meta)
    pklfile=osp.join(self.outdir,filename)
    pd.to_pickle(pklfile)
    return
