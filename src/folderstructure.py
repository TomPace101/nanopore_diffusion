"""Define the expected folder structure

This module defines a dictionary ``folder_structure`` used by some classes in the ``locators`` module to determine file locations."""

#Standard library
import copy
import os
import os.path as osp
import sys

# #Locate source folder
# if 'SRCLOC' in os.environ.keys():
#   srcfolder=osp.normpath(osp.abspath(os.environ['SRCLOC']))
# else:
#   srcfolder=osp.abspath(osp.split(__file__)[0])
# 
# #Add the source folder to the python path if not already present,
# #to allow importing modules
# #add python code folder(s) to path
# if not srcfolder in sys.path:
#   sys.path.append(srcfolder)

from filepath import Path

# #Use filepath.Path now
# srcfolder=Path(srcfolder,isFile=False)
srcfolder=Path(__file__).parent

# simulator_modules_folder=srcfolder / 'simulators'
# if not simulator_modules_folder in sys.path:
#   sys.path.append(simulator_modules_folder)

#Locate data folder
if 'DATALOC' in os.environ.keys():
  # datafolder=Path(osp.normpath(osp.abspath(os.environ['DATALOC'])))
  datafolder=Path(os.environ['DATALOC']).expanduser().reslolve()
else:
  datafolder=srcfolder.parent / 'data'

#Default folder structure
default_folder_structure={
  'CustomizationFile':['customizations'],
  'RequestFile':['requests'],
  'MeshGeomdefFile':['mesh','geomdef'],
  'MeshTemplateFile':['mesh','templates'],
  'MeshGeoFile':[0,'mesh',1],
  'MeshMshFile':[0,'mesh',1],
  'MeshXmlFile':[0,'mesh',1],
  'MeshHdf5File':[0,'mesh',1],
  'MeshGmshOutFile':[0,'mesh',1],
  'MeshMetadataFile':[0,'mesh',1],
  'SolutionFile':[0,1],
  'PostprocFile':[0,1]
}

#run-time alterable folder structure
folder_structure=copy.deepcopy(default_folder_structure)
