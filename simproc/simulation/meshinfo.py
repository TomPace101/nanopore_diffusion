"""FEniCS mesh support"""

#Standard library
import os.path as osp

#Site packages
import fenics as fem

#This package
from ..requesthandler import yaml_manager

class MeshInfo:
  """Bunch of mesh-related data

  Attributes:

    - mesh = FEniCS Mesh
    - facets = FEniCS MeshFunction of gmsh Physical Surface number (3D) or Physical Line number (2D)
    - cells = FEniCS MeshFunction of gmsh Physical Volume number (3D) or Physical Surface number (2D)
    - metadata = dictionary of metadata about the mesh, such as parametric locations

  A note on the terminology used in FEniCS and gmsh:

  |  The FEniCS information below is from page 185-186 of the FEniCS book.
  |  d = number of dimensions in entity,
  |  D = number of dimensions in problem (maximum entity dimension)
  |  D-d = "codimension" of entity
  |  Terms:
  |    D=2, d=1: fenics facet (facet_region xml) = fenics edge = gmsh physical line
  |    D=2, d=2: fenics cell (physical_region xml) = fenics face = gmsh physical surface
  |    D=3, d=2: fenics facet (facet_region xml) = fenics face = gmsh physical surface
  |    D=3, d=3: fenics cell (physical_region xml) = fenics ____ = gmsh physical volume
  |    also, d=0 is a fenics vertex"""

  def __init__(self,mesh=None,facets=None,cells=None,metadata=None):
    if mesh is None:
      self.mesh=fem.Mesh()
    else:
      self.mesh=mesh
    if facets is None:
      self.facets=fem.MeshFunction("size_t", self.mesh, mesh.geometry().dim()-1) #Facets are of dimension d-1
    else:
      self.facets=facets
    if cells is None:
      self.cells=fem.MeshFunction("size_t", self.mesh, mesh.geometry().dim()) #Cells are of dimension d
    else:
      self.cells=cells
    if metadata is None:
      self.metadata={}
    else:
      self.metadata=metadata

  @classmethod
  def load(cls,mesh_hdf5,meshmetafile,loadfuncs=True):
    """Load Mesh and MeshFunctions from HDF5 file, and mesh metadata from yaml
    
    Arguments:
    
      - mesh_hdf5 = path to the hdf5 file for the mesh
      - meshmetafile = path to the mesh metadata yaml file
      - loadfuncs = optional, True (default) to load mesh functions, False otherwise"""
    #Load mesh metadata file, if it exists
    if meshmetafile is None:
      metadata=None
    else:
      metadata=yaml_manager.readfile(str(meshmetafile))
    #Initialize empty mesh
    mesh=fem.Mesh()
    #Open HDF5 file
    hdf5=fem.HDF5File(mesh.mpi_comm(),str(mesh_hdf5),'r')
    #Get the mesh
    hdf5.read(mesh,'mesh',False)
    #Initialize the object
    self=cls(mesh=mesh,metadata=metadata)
    #Load meshfunctions if requested
    if loadfuncs is True:
      hdf5.read(self.facets,'facets')
      hdf5.read(self.cells,'cells')
    #Done
    hdf5.close()
    return self
