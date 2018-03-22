
"""Check that the required software is available."""

import subprocess
import sys
assert sys.version_info.major == 3, "Python 3 required."

#check required python modules
import yaml
import jinja2
import numpy
import scipy
import matplotlib
import pandas
import fenics as fem

#fenics version check
target_fenics_version='2017.1.0'
fenics_version_msg_template="This code was written for FEniCS version '%s'. Detected version '%s'."
assert fem.DOLFIN_VERSION_STRING == target_fenics_version, fenics_version_msg_template%(target_fenics_version,fem.DOLFIN_VERSION_STRING)

#check that gmsh runs
gmsh_cmd="gmsh --version"
retcode=subprocess.call(gmsh_cmd,stderr=subprocess.DEVNULL,shell=True)
assert retcode==0, "gmsh failed to execute correctly"

#tests passed
print("OK")