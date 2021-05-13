
- contour plotting: see notes 2020-10-29, 2020-10-30
- clean up rst formatting in module docstrings (see warnings from the sphinx build)
- example of what the yaml looks like for a request with no needed arguments. Clean logs is a good example. But do mention that a name is usually a good idea.
- Log file cleanup requests require their own log file to be set up in order to find the logging directory. (a file which, of course, will be auto-erased).
- check the module dependency file: in particular, I'm not sure logging is listed everywhere it should be
- MPI request children don't get logged: logging parameters not passed? (2020-09-28)

singularity recipes:
- add sphinx to make a "complete" recipe, and take testing out of the "minimal" one (`fenics_2019`).
- add git

_FEATURE_ collection: a way to add data from the job list columns into the new table
For example, add an attribute: job_columns,
which is a mapping from column names in the job table to the column name in the new table.
This would have made things easier for me several times already.

_FEATURE_ better python api
This is really the biggest problem with simproc.
It's set up assuming you'll be accessing it from yaml.
There are lots of useful functions that should return values and accept arguments.
But instead of being able to do it that way,
I have to store data into attributes and get the results from some other attribute.
The python api should make these available as functions.
The yaml api should then make use of that python api.
I think that perhaps one particularly bad consequence of this
is when things like parallel task execution
require writing out separate yaml files.
Ultimately, a request is a function call.
It's just a call where the arguments and return are stored,
so I can check to see if I need to actually call it again or not.
It's a big memoization, basically.

_FEATURE_ replace/improve doit
(This is a new summary of `Dependency tracking` far below.)
As much as a I like doit, it does have some issues:
- the signatures are stored in a file that isn't human readable
- there is no way to manually force a signature update
- signatures aren't updated on a task-by-task basis, so if doit crashes before all tasks finish, all the intermediate progress is lost
- tasks can't modify the task list, complicating auto-generation of the task list in nontrivial circumstances
- there's not an easy solution to this one, but: when my task names change, the database loses any previous state and assumes the task is new, because the task names are the database keys. It would take pattern matching to fix this, I guess?

**note that some items below have already been resolved**

**note that items related to generated requests are probably made obsolete by the job list module**

_ISSUE_ validation: debug.yaml
The GeneratedVariationsRequest prevents the validation from being run with doit.
Appending to the same file means two different requests have the same output file.
We need to revise this test so the generated variations have different output filenames.

_ISSUE_ run generated requests simultaneously
Does this work now?
(Potentially no longer relevant, based on changes to request generation: see below)
Once you get it working, include an example in the validation.
Maybe this is a secondary class,
which creates the queue for the existing class,
as the child of the secondary class.
Related issue: some requests really do need to be handled sequentially.
Or, more generally, you need to be aware of dependencies when you try to run in parallel:
Some things have to finish before others can start.

_ISSUE_ with new request generation approach
The new approach is based on using templates of yaml files,
like the way `paramgen` used to work.
But it actually still doesn't work as well as that did.
It's still a two-step process: generate the file, then run it.
Another, minor issue: no way to generate sequential ID for the sub-requests.

_ISSUE_ task dependencies
doit got it wrong pretty bad in redux_electrolyte once
The plot was wrong because things hadn't been rerun.
How did this happen?
It was when I was trying to do the second curve.
I reverted, and somehow it didn't re-run the smoluchowski solutions.
It did re-run the extraction, but it must have missed the change in the configuration
of the smoluchowski solutions somehow.
More generally, most of the validation doesn't test the accuracy of the input and output files lists.
Or the configuration information.
How can we check that this is correct for all the different request types we have?

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
- generating MPIRun requests

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

_ISSUE_ should the stuff in locators that allows control from yaml be moved to another module?
Those should probably have schemas, right?
Then they need to load that module, but schema needs the locator classes for validation.
Actually, maybe UpdateFolderStructure doesn't need one.
But the others do, it seems.
Still, it just seems like configuring/setting locators from yaml
is separate from the base functionality of locators themselves.

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
That should be working now.
Except that the last time I tried, the MPI simulation crashed without getting a solution!
(See log 2019-07-18.md)

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

_FEATURE_ a variant of parallel request that does one item first, then does the rest in parallel
This is to help avoid FFC cache collisions.
I have now attempted this in simultaneous.py, but have not tested to see if it actually works.

_ISSUE_ Templates for request generation may need to contain child requests
The reason the templates are dictionaries (instead of just request instances themselves)
is that a template may not define enough information for a request to be valid by itself.
It may need to be filled in during the generation process for that to work.
But that logic would apply to children of the template as well.
Currently, that won't work.
Children of a template have to be provided as request instances in the yaml file.

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
First attempt now implemented as nested.Stored
I still can't help but wonder if this is what the "descriptor" protocol in python was made for.
https://docs.python.org/3/howto/descriptor.html

_CONCERN_ in `collection.py`, `CollectionRequest` calls the `complete` method of `RawCollectionRequest`
Maybe it should actually create a single child instead.
This is marked as a TODO item in the file.

_FEATURE_ collection requires all rows in the table to have data in each column
We could remove this requirement if we identified a key column, which then had to appear in each mapping.
Of course, your simulations have to use the same values for this key, and output them to the results files.

_ISSUE_ `CollectionRequest` and `SimpleCollectionRequest` don't seem useful in practice.
I don't always know how to get handles to all the requests they would need.
As of now, I've been explicitly listing all the files, which is a pain.
Maybe if I could list the request names, and let them be retrieved from the named object store,
that would be different.

_ISSUE_ in collection, there needs to be some way to simplify the mapping specification
Right now, if all my entries in the file are suitable column names,
it doesn't matter.
I have to list them all anyway, as themselves!
(e.g. `column_1: column_1`)

_FEATURE_ see TODO items in `generate.py`

_FEATURE_ `simrequest.py` had a workaround for the lack of the `schema` module for checking the "conditions" attribute.
Update this to use the schema module instead if possible.
And, of course, the other simulator modules follow this same approach.

_FEATURE_ request generation for a list of requests, instead of just one.
Allow a list of requests instead of just one.
Basically like having a for loop.
For that, you need loop variables that can be used within the loop.
This could potentially be done with attribute paths (now the Stored class),
but only if children can find their parent.
Right now, they can't.
For that matter, maybe it's not just parents that need to be found.
It could be any related request from which I might need to copy some attribute.
So, we want a dictionary of related requests.

Related: why specify a request_type, instead of just making the template itself be of that type?
The answer is that the template may not have enough information to be a fully valid request on its own.
It might need the information put into it from these other sources for that.

Related: finding files located by other requests.
Except that the other request is the one that needs to render the locator.

See logs:
  - 2019-08-05.md
  - 2019-10-07.md

The checkerboard problem is now a pretty good example of how to do this manually.
It runs a bunch of requests, with a parent request to define them.
Maybe some things from that could be generalized?

# Dependency tracking
Here are some limitations of doit that I'd like to resolve:
- If I should get something to run without doit, I have to rerun it to get its status updated.
  There should be a way to manually use the current state of a task and tell it "this is up to date as of now".
- Sometimes I have different sets of "old" files for meshes, solutions, and post-processing.
  I can't tell which one went with what. 
  It would be nice if the file hashes were stored in a human-readable format for inspection.
  (Both input and output file hashes for each task, updated on task completion.)
  And not just the current versions; a history of them might be helpful as well.
  Then I could identify old stuff.
  Some of this could also be potentially worked around by using git-LFS to "track" the output files.
  Then I could see their versioning.
- If a task crashes so hard it takes doit with it, none of the task status updates are written to disk.
  I have had this happen with tasks that produce memory allocation errors.

Maybe the way to do this is to implement some things in my own fork,
then submit pull requests.
But to make that work, I'd need to switch to using my fork inside
the singularity images.
