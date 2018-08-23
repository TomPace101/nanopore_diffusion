#Import all modules that define classes loadable from yaml
from . import filepath
from . import locators
from . import requestfile
from . import customization
from . import debug

#Modules for doctests
doctest_modules=[filepath, debug]
#Files for doctests
doctest_files=[]
