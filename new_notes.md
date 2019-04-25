
_ISSUE_ gmsh version change has broken some mesh file comparisons.
This was already an issue with the hdf5.
We need a better way to do this.
Maybe have a file size comparison list, like we have a list for file comparisons.
Maybe even set up deeper inspection classes for particular file types:
xml, yaml, hdf5.
You'll have to run on another computer to find the issues.

_TODO_ test out line_profile

_ISSUE_ There are input files for the simulation test that aren't tracked in git.
We need to have another validation step generate this file first.
That will require the expression projector.

_ISSUE_ the simulation test doesn't compare the file output
Maybe we even need a way to compare yaml files other than just bytewise.
For example, compare floats to a limited precision.

_ISSUE_ the homogenization simulation requests have no property schema, and an incomplete docstring

_TODO_ in writefield_outputfiles, detect MPI and list pvtu files

_ISSUE_ we really should combine the 2D and 3D homogenization files once we get them working, to reduce duplicated code
The PeriodicBoundary classes will have to stay separate.
This requires a way to identify the dimensionality of a mesh.
We could put it in the mesh metadata.
Use `mesh.geometry().dim()` to get number of dimensions in a mesh.

_ISSUE_ the mesh locators don't match the case of the others
For that matter, are there names consistent with corresponding attributes elsewhere?

_ISSUE_ Simultaneous requests may create a new directory for the temporary request input files.
But cleanup won't remove this directory if that's the case.
Maybe requests need to be aware of their own temporary files as well?

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
  - DONE Basic simulation as subclass of customizable request, with the key methods.
  - DONE MeshInfo (similar to what's in simulator_general now, but not taking modelparams)
  - General input tables? (see `p20180819_InputTable`, or should use pandas?)
  - Specific input tables for species, domains, and species-in-domain
  - DONE Weak Form support as exists in simulator_general now
  - Library of common weak forms
  - expression projector?
  - LPB/Smol simulator
  - unhomog fickian diffusion?
- Post-processing:
  - Collection request to generate table, with customization to control how simulation output maps to columns, and maybe specific data for select fields as well
  - Plot request similar to what's already there, but with customization as well
  - zipping/unzipping: see below
- Validation:
  - DONE Validate that output files are as expected (new request type, which is tested by doctest)
- Request generation:
  - DONE Requests store themselves in a yaml file
  - DONE Request that can parametrically generate child requests
- Customization:
  - DONE a request that can monkey-patch itself
  - allow user to specify python files containing classes that can be added to yaml registry

# New Features/Improvements

_DESIGN ISSUE_ Paths: you never know if it will be a locator, path, or string.

One idea was to try to switch all locator rendering to the late method.
The challenge is that sometimes you're working with locators that belong to a different request.
So you've got to watch out for that, and have access to the request, not just its locators.
(This is actually closer to the way it worked before: ParameterSet has methods for constructing file paths.)
Maybe it won't work, but you can try.

But does that really resolve the issue?
Render-at-point-of-use will work for the locator-to-path step.
But still I don't know if I have a path or a string.
I can assume it's a string, but then I can never use paths.
I can cast it to a path, but that might be redundant.
Right now, there are possible cases where the code requires one or the other,
but there isn't a conversion done.
How do we systematically find such cases?
You could do this:
- render to make sure you have a path or string
- then case to path or string, as needed.
But that's tedious.
And prone to being forgotten.

_FEATURE_ run with doit without dodo.
See old notes about this.
This requires digging into doit and copying out some of its code.
The advantage is obvious: less time waiting for things to rerun.

_FEATURE_ should attribute paths be moved up to request itself?
(This is `get_nested` and `set_nested` in simrequest.py)
It could be useful in collection, too, actually, so maybe it should be its own module.
This gets back to the idea of "memory locators".
Maybe we need better names than get_nested and set_nested.
Maybe we need a better name than "attrpath" ("attribute path").
For one thing, that makes it sound like its a filesystem path.
At the very least, we need to explain this somewhere.
Maybe even have an example in the tutorial.

_FEATURE_ a variant of parallel request that does one item first, then does the rest in parallel
This is to help avoid FFC cache collisions.

_FEATURE_ Why aren't shell requests customizable?

_FEATURE_ confirm input files exist as part of pre-run check
I wrote code to do this, then I realized that not all the input files are actually required.
How can the code make that distinction?
Maybe it can't.

_FEATURE_ it would be better if meshinfo could query the HDF5 file about its components
rather than requiring a keyword argument.
