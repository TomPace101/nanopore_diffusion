#Import all modules that define classes loadable from yaml
from . import cleanup
from . import comparison
from . import customization
from . import debug
from . import filepath
from . import generate
from . import locators
from . import mpi_run
from . import nested
from . import request
from . import requestfile
from . import shell
from . import simultaneous
from . import templates
from . import yaml_manager

#Modules for doctests
doctest_modules=[filepath, debug, comparison]
#Files for doctests
doctest_files=[]
