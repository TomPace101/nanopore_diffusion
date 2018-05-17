"""Convert XML mesh and mesh functions to HDF5 format"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp
import sys

#Site packages
import fenics as fem

#Local
import folderstructure as FS
import common

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

##TODO!!
