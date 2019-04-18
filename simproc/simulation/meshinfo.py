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

  def __init__(self,mesh_hdf5,meshmetafile,loadfuncs=True):
    """Load Mesh and MeshFunctions from HDF5 file, and mesh metadata from yaml
    
    Arguments:
    
      - mesh_hdf5 = path to the hdf5 file for the mesh
      - meshmetafile = path to the mesh metadata yaml file
      - loadfuncs = optional, True (default) to load mesh functions, False otherwise"""
    #Initialize empty objects
    self.mesh=fem.Mesh()
    self.facets=fem.MeshFunction("size_t", self.mesh)
    self.cells=fem.MeshFunction("size_t", self.mesh)
    #Read in data from HDF5 file
    hdf5=fem.HDF5File(self.mesh.mpi_comm(),str(mesh_hdf5),'r')
    hdf5.read(self.mesh,'mesh',False)
    if loadfuncs is True:
      hdf5.read(self.facets,'facets')
      hdf5.read(self.cells,'cells')
    hdf5.close()
    #Load mesh metadata file, if it exists
    if meshmetafile is not None:
      self.metadata=yaml_manager.readfile(meshmetafile)
