
_ISSUE_ locators don't have descriptions of their purpose.
This is important: people won't know how to use them without this.
There should also be a good way to see all defined locators.
(You can see them all with `locators.folder_structure`, but I mean something more helpful, somehow.)

_ISSUE_ gmsh version change has broken some mesh file comparisons.
Mesh xml files have variations in the ordering of indices, and lower decimals of floating points
There was already a similar issue with the hdf5.
We need a better way to do this.
Maybe have a file size comparison list, like we have a list for file comparisons.
Maybe even set up deeper inspection classes for particular file types:
xml, yaml, hdf5.

_TODO_ missing validation:
- line_profile
- collection
- defining a series from a collected data file (doing this will also test CommandSequenceRequest, indirectly)
- plotting

_ISSUE_ There are input files for the simulation test that aren't tracked in git.
We need to have another validation step generate this file first.
That will require the expression projector.

_ISSUE_ the simulation test doesn't compare the file output
Ideally, these simulations should be based on problems where we can compare to an analytical solution.
Maybe we even need a way to compare yaml files other than just bytewise.
For example, compare floats to a limited precision, or other form of tolerance.

_TODO_ in writefield_outputfiles, detect MPI and list pvtu files
There are also vtu files generated that don't match the current naming expectation in this case.
The workaround is to not write out pvd files in MPI.

_ISSUE_ we really should combine the 2D and 3D homogenization files once we get them working, to reduce duplicated code
The PeriodicBoundary classes will have to stay separate.
This requires a way to identify the dimensionality of a mesh.
We could put it in the mesh metadata.
Use `mesh.geometry().dim()` to get number of dimensions in a mesh.

_ISSUE_ Simultaneous requests may create a new directory for the temporary request input files.
But cleanup won't remove this directory if that's the case.
Maybe requests need to be aware of their own temporary files as well?

_ISSUE_ It's been a while since we built documentation with sphinx.
We probably want to start from scratch on that for the refactoring.

_ISSUE_ when requests are written to yaml (MPIRunRequest and SimultaneousRequestQueue),
could there be modules that need to be loaded by the new process as well?
Could there be folders that need to be added to the python path?
You could have the customization module track this.
When its relevant classes are used, it keeps track of the additions in a module variable,
and then can write those items back out to yaml.

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
  - DONE: LPB/Smol simulator
  - unhomog fickian diffusion?
- Post-processing:
  - DONE: Collection request to generate table
  - DONE: Plot request
  - zipping/unzipping: see below
- Validation:
  - DONE Validate that output files are as expected (new request type, which is tested by doctest)
- Request generation:
  - DONE Requests store themselves in a yaml file
  - DONE Request that can parametrically generate child requests
- Customization:
  - DONE a request that can monkey-patch itself
  - DONE allow user to specify python files containing classes that can be added to yaml registry (this means loading the module)
  - DONE allow user to specify folders to be added to the python path, so other python files can import their modules without the data files listing them

# New Features/Improvements

_FEATURE_ Run simulations with MPI
The issue is data extraction: each process only has part of the mesh.
(See log 2018-05-29.md)

fenicstools
https://github.com/mikaem/fenicstools/wiki

So, what will it take to get this working:
- rebuild singularity images with fenicstools installed
- use fenicstools Probes to extract data
- MPI gather results from different processes (or does Probes do this for you?)
- only rank 0 process should write output file

For now, we're working around this by doing extraction in single-process mode.

_FEATURE_ run with doit without dodo.
This requires digging into doit and copying out some of its code.
The advantage is obvious: less time waiting for things to rerun.

You can call doit from within python; see doit.run
http://pydoit.org/cmd_run.html#using-the-api
Probably you'll want to examine that code,
and do something similar.
You can't use it directly, because it calls `sys.exit`.

What format to use this?
Could add a command line option to "run" with doit.
But sometimes I just want to see task states.
We don't want to replicate all the doit command line stuff.
The dodo file does work for that.
Instead, the goal is to be able to do those things from within python.

(It would be nice if doit had a better python api by itself.)

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

_MAYBE_ the mesh locators don't match the case of the others
For that matter, are their names consistent with corresponding attributes elsewhere?
Same thing for "modulefile" in `customizations.py`.

_ISSUE_ Unconfirmed: once a child request of an MPIRunRequest crashed,
but it seemed like the parent request kept waiting.
But I can't get that to happen again.

_FEATURE_ attribute/item "locators"
Needs a better name to distinguish it from file path locators,
but the basic idea is that anywhere in my code that currently accepts
a "nested attribute path", should really be able to accept the data directly instead.
But, we then have a class that works similar to a locator:
it tells you where to find the data.
Just like locators, there's a method to "render" objects,
which just returns the direct data if that's what was provided.
In fact, the objects themselves don't need a method for this;
it's all up to the owning object to handle properly.
But then what about where to "store" the result?
That's not really a use case for this.
When you need to store a result, you know you have to tell where to put it.
The use case is when you get a value,
you can either provide the value directly,
or provide the location of the value.
Terminology:
"owner": the object needing the value (which may also possess the value as a nested attribute)
"attribute locator": the object that tells the owner where to find the needed data
"direct value": not an attribute locator; for use by the owner without alteration
"argument": something the owner gets passed that could be either a direct value, or an attribute locator
"characteristic": something that distinguishes attribute locators from direct values.
How it should work:
The owner checks if the argument possesses the characteristic.
If it does, it calls owner.get_nested.
If not, the value is returned directly.
This could all happen, for example, in owner.render.
All we need to do, then, is decide what the characteristic is.
Maybe it's just "isinstance"?
The place to put the attribute locator would be in nested.py.

So now I made a first pass at an implementation of this.
But I'm not sure I like how it works.
I couldn't use "render", because that has to be guaranteed to return a file path,
not an arbitrary value.
Also, what would you do if you had a locator stored somewhere?
You want to retrieve it without rendering it.

Nothing else uses it yet, so it's dead code at the moment.

I still can't help but wonder if this is what the "descriptor" protocol in python was made for.
https://docs.python.org/3/howto/descriptor.html

