
_TODO_ go ahead and switch request refactoring into the master
First, make a branch tracking what is now the master branch.

_TODO_ There are input files for the simulation test that aren't tracked in git.
We need to have another validation step generate this file first.

_FEATURE_ a variant of parallel request that does one item first, then does the rest in parallel
This is to help avoid FFC cache collisions.

_FEATURE_ Why aren't shell requests customizable?

_FEATURE_ confirm input files exist as part of pre-run check
I wrote code to do this, then I realized that not all the input files are actually required.
How can the code make that distinction?
Maybe it can't.

_ISSUE_ There's something weird in `request.py`.
`_compile_file_list` uses attribute `taskname` which doesn't exist, but I've never gotten an error about this.
`task_definition` says it will not return a task if the taskname is None, but doesn't do any such check.
`assure_output_dirs` doesn't use the outputfiles property, because that property is supposed to return those of children as well.
Now `confirm_inputfiles` does the same thing.
Maybe we need a way to specify input/output files for this request alone, or this request with children.

# Refactoring

- DONE Top handler: a new module that serves to dispatch requests to the handlers defined in other modules
- Base classes:
  - DONE Request - abstract only: defines the interface
  - DONE Request to run all of the requests in another file. Or maybe even a way to specify a subset of them. (see note about doing this on command line)
  - DONE Request to run a script (or any shell command)
  - DONE Template requests: requests that produce a file from a template and input data
  - DONE Cleanup requests: requests to delete the output from previous runs of requests
  - DONE Customizable requests
  - DONE Request to run child requests in parallel (that is, simultaneously)
  - DONE Request to run a child request with mpirun
- Data Locators:
  - DONE DataFile: defines a file relative to the data folder
  - DONE update folder structure from python code or yaml file
  - DONE update datafolder location from code or yaml file
  - ONGOING others are to be created along with the requests they relate to
- Mesh generation:
  In the info below, reqdata means the request includes data other than just the input and output files
  - construct template from geometry definition (reqdata -> .geo.jinja2) [not the way it was done before. will this work?]
  - create geo file from template and values (.geo.jinja2 + reqdata -> .geo) [OR is this too case-specific? The request data->template input conversion step varies.]
  - DONE run gmsh (.geo -> .msh)
  - DONE run dolfin-convert (.msh -> .xml (3))
  - DONE hdf5 conversion (.xml (3) -> .hdf5) [write a python function, in a module supporting command line]
  - DONE later steps only (.geo -> .hdf5)
  - all but template construction (.geo.jinja2 + .yaml -> .hdf5)
  - all steps (.yaml -> .hdf5)
- Simulation:
  - Basic simulation as subclass of customizable request, with the key methods.
  - MeshInfo (similar to what's in simulator_general now, but not taking modelparams)
  - General input tables? (see `p20180819_InputTable`, or should use pandas?)
  - Specific input tables for species, domains, and species-in-domain
  - Weak Form support as exists in simulator_general now
  - Library of common weak forms
  - MORE
- Post-processing:
  - Collection request to generate table, with customization to control how simulation output maps to columns, and maybe specific data for select fields as well
  - Plot request similar to what's already there, but with customization as well
  - zipping/unzipping: see below
- Validation:
  - Validate that output files are as expected (new request type, which is tested by doctest)
- Request generation:
  - UNTESTED Requests store themselves in a yaml file: should work already, just test it
  - DONE Request that can parametrically generate child requests
- Customization:
  - DONE a request that can monkey-patch itself
  - allow user to specify python files containing classes that can be added to yaml registry

# BUGS
- Simultaneous requests may create a new directory for the temporary request input files.
  But cleanup won't remove this directory if that's the case.
  