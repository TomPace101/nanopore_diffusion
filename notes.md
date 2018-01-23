

_TODO_ create a path object (unless pathlib can do it) supporting:
- folder
- base name
- extension
- file name (property)
- full path (property)
This would simplify things a lot if we had it.
Or would it?
You have to use this class as the appropriate object attribute,
meaning at initialization it has to be there.
That doesn't fit with the way ParameterSet works.

_TODO_ in buildgeom, validate geometric inputs (different formulas for different geometries)

# Refactoring

- let objects provide their own info needed by doit tasks in a standard structure
  - DONE: changes within useful
  - Load-level Attributes/properties/methods needed (see below)
    - \_child_attrs (class attribute, list, if needed)
    - \_inputfile_attrs (class attribute, list, if needed)
    - \_more_inputfiles (class attribute, list, if needed)
    - \_outputfile_attrs (class attribute, list, if needed)
    - \_more_outputfiles (class attribute, list, if needed)
    - \_folders (class attribute OR instance attribute, dictionary)
    - \_required_attrs (class attribute, list)
    - \_config_attrs (class attribute, list)
    - \_taskname_src_attr (class attribute), string
    - run()
    - and the __main__ portion of the module as well (change process_function to process_method, usually omitted)
    - change 'from folderstructure import \*'
  - Objects needing each of these:
    - DONE: paramgen.ParameterGenerator
    - DONE: buildgeom.MeshParameters
    - new module to generate .msh files from .geo, based on yaml (run_gmsh)
    - new module to generate .xml files from .msh, based on yaml (run_dolfin_convert)
    - solver_general.ModelParameters
    - postproc.PostProcParameters
- modify dodo.py and the various task_ files to use the new approach, reusing (after extracting out) code from the command-line version where possible
- document the general approach in README.md

_TODO_ shell command files for generating .msh and .xml files, which doit then calls?

_ISSUE_ doit tasks can't generate meshes when there are no models
This was kind of as intended: it only generates meshes that are needed.
But maybe you should generate all the meshes for the basenames requested in control.yaml.

Currently, control.yaml actually lists the model files to be read,
not the mesh files.
That way, the mesh files and model files are not required to have the same names.

If you want to change that, you need to think carefully through all the consequences.
The easiest way to deal with the issue here is probably to have a "null" model of some kind,
that will not generate any tasks.
When you only want to generate meshes, use that model.

Or, we could do it this way:
control.yaml contains a list of either model or mesh (and maybe even postproc or paramgen) files.
You run the ones that exist, and don't complain about the ones that don't.
(Unless no files are found; then you should issue a warning.)
That seems like it could cause errors: a typo in a filename wouldn't be as detectable.

_FEATURE_ doit tasks for parameter generation
But other tasks are generated based on reading the output of parameter generation tasks.

doit has a way to resolve this, of course:
http://pydoit.org/task_creation.html#delayed-task-creation

That doesn't completely resolve the issue:
the code that generates the other tasks is outside any task function.
So even if the other tasks happen after parameter generation,
the parameters for them have already been read in.

Perhaps there could be a task to populate those global parameters
(that is, outside tasks they initialize empty, and tasks fill them in),
which is delayed until after parameter generation.
The other tasks are delayed until after population.

Or maybe it doesn't need to be a task at all.
Maybe you just include the code to generate the parameters as 'always execute'
in the dodo file, then generate the tasks.
This is somewhat wasteful because parameter generation could take time.
And it would cause file dates to reflect last run,
not necessarily the last time something actually changed.
(Though you shouldn't always count on that anyway.)

__ISSUE__ Some of the current problems suggest we are too tightly coupled to doit.
- The changes needed for 2D geometry
- Not being able to generate meshes for which there are no models
- difficulties with parameter generation tasks
- not having a way to easily generate .msh and .xml files without doit
- challenges in using MPI for fenics but not everything else

What should be happening here is just that doit runs things a human being could conceivably do,
but in an automated and organized fashion.
That is, it just ties together command-line programs.
Of course, if they are python modules/scripts, it can call functions inside them instead.

Or maybe that's the problem.
In trying to set things up so doit can use the objects directly,
I've set up two parallel systems of doing things.
And sometimes they conflict.

Take a look at the post-processing tasks, for example.
You couldn't possibly run them without doit to figure out what to call.

Currently:
- paramgen can be called directly to generate yaml input files from yaml
- buildgeom can be called directly to generate .geo files from yaml
- individual solver modules can be called directly to solve models from yaml

A more consistent, comprehensive approach:
- paramgen as now
- buildgeom as now
- something to generate .msh files from .geo, based on yaml
- something to generate .xml files from .msh, based on yaml
- something to find and run the appropriate solver, based on yaml
- something to run post-processing tasks in a yaml file, which also specifies the folder where the models can be found (this is a new parameter for that file)

The catch is that we lose some granularity over re-running things.
Right now, we can re-run individual models from a long list.
The new approach would re-run the entire list if it is specified.
That's how we got to the current approach.
The tight coupling is needed to figure out what really needs to be run and what doesn't.

Is there a way to get both?
The command line programs would need to be able to take additional arguments
that specify wich documents in a multi-doc yaml file should be run.
They need to know a property by which to identify such documents.
So an argument for the property name, and another for its value.
These are optional, but both must be provided.
(In argparse this is done by an option that takes two arguments:
  https://docs.python.org/3/library/argparse.html#nargs
  e.g. '--select modelname debug001'
  or you could let it take at least two arguments, and all beyond the second are additional cases to run)

How does that help?
It sets up a more consistent entry point between doit and the command line.

And so, really, the goal here is to have less duplicated code between the two.
That is, there should only be one place where the objects are loaded from yaml,
not one for doit and a different one for the command line.
But the way we determine file and parameter dependencies is by loading the objects
and examining them.
Is there a way to do it without that?
For example, maybe doit could just leave things as dictionaries, or namespaces.
Yes, that's how to somewhat decouple without losing the granularity.

Also, the doit task generation should be more consistent between the different tasks.
What does that look like?
- for each stem, see if there is a file appropriate for this task
- if so, load all the documents
- do the investigation to see which ones are not up-to-date
- generate tasks, then with the entry point as the python action,
- OR a command action with only required tasks specified by the --select argument

For that matter, command line execution should be more consistent as well.
- accept one required yaml file, and optionally the select argument above
- then load the yaml file as a sequence of namespaces, extract the selected documents (or all docs if no selection)
- call the task-specific entry point with the list of namespaces
- the entry point converts namespaces to objects
- delete the list of namespaces

Actually, we should just use the native objects from loading the yaml (ie dictionaries) instead of namespaces.
And should it be that the entry point does only a single task?
The command line calls the entry point once for each document to process, with the dictionary as the only argument.
Doit:
- loops over the steps (the modules that are runnable from the command line), and in each step
- loops over the filenames in control.yaml, and for each filename
- sees if there is such a file for this step, and if so
- it loads all the documents,
- sees which ones are out of date, and for each of those
- creates tasks which call the entry point with the document dictionary

And you know, we could actually still use the appropriate objects,
if for meshparameters we just created a new parent attribute for the geometric parameter values.
That just has some consequences of its own.
- parameter generation changes (fairly minor)
- post-processing changes (somewhat extensive, but easy to find)

Yes, this means two separate loadings of the document.
But we taught these objects to load themselves from yaml for a reason.

Now we get to the real problem: loading meshparameters for each model.
It's inefficient to do it one model at a time, reloading all the meshes in the specified file each time.
It's slow enough already.
That's what got us into loading all the meshes needed by all the models first,
and storing them in a dictionary.
(This is a time-memory tradeoff.)

One way to do this is to have the entry point function accept not just one object,
but the generator itself. As originally planned.
But this is where the conflict with doit happens.
Each object needs its own task.

Another way is to memoize (or, really, just cache) the results of the mesh parameter file loading.

So here's how duplicate code between the command-line and doit versions happens:
They each need the same information (input file locations, parameters, output files)
but they do different things with them.
So we need to create little functions to generate these things.
That way, both doit and the command line can call the same functions,
rather than having to duplicate the same code.
Maybe these should even be methods of the objects.
For example, the PlotFigure objects have locate_data methods.

What do they need to provide:
- a dictionary of their file dependencies, with useful keys: inputfiles
- a dictionary of their targets, with useful keys: outputfiles
- a name that could be used for a task name (maybe we already have this, just use the modelname, meshname, etc. But we need a way to say what this is. So it needs to be a property.)
- a method that can be used as an action: run()
- a method to return a configuration string suitable as an argument to config_changed: config()

The "run" method might obviate the need for some of the functions I created today.
That is, that code might be better as a method in some cases.

Maybe some of this functionality could be implemented in useful.ParamterSet.
That is, give it slots for inputfiles and outputfiles, and initialize them as dictionaries.
Maybe they should be properties, so that they can be computed when needed?
Or should they be loaded at initialization?

So a typical doit task definition would look like:
{'name': self.taskname,
 'file_dep': list(self.inputfiles.values()),
 'uptdodate': config_changed(self.config,
 'targets': list(self.outputfiles.values()),
 'actions': [(self.run,)]}

In fact, maybe the objects can just return their own doit task definition.
Except that meshparameters needs 3 tasks, not 1.
Or, we create new objects and use those, but they are read from the same yaml file.
Or, the method that returns task definitions is a generator.
In the simple case, it only yields one task definition.
But for MeshParameters, it can yield 3 of them.
No, that doesn't work for the command line.
The action's wouldn't all just be "run".
For the typical command line approach, we need new modules, with their own classes.
Still, it seems less efficient.
We have to load the same yaml file 3 times, instead of only once.
Maybe that is an acceptable price to pay for consistency?
Maybe instead of run, the object gives a list of possible actions.
You select the one you want from the command line,
and optionally can omit this to run all of them.
So now the task definition can be a generator after all.
useful.run_cmd_line needs to accept an optional command line argument
specifiying a task to be done with the loaded object.
So the objects have a dictionary {cmd_line_task_specifier: method_name}.
This dictionary is called "actions".
But that complicates the inputfiles and outputfiles as well.
(Different actions use different ones.)
Instead, we do need the 3 different objects.
So, there is no reasonable way around loading the file 3 times,
as each object will be loaded from the file.

Also, .geo files aren't built from one object, but two:
- MeshParameters
- GeometryDefinition
So maybe there needs to be a single object that has both of them as attributes?
Or maybe MeshParameters needs to allow a GeometryDefinition as a child attribute.

Useful.run_cmd_line then would no longer need a function to call.
It could just call the run method of the object it is passed.

Except that wouldn't quite work for buildgeom.
Unless MeshParameter WAS the object that was run,
and GeometryDefinition was just a child of it.
Of course it should be.
And the functions should be methods, usually of MeshParameters.

The memo stuff in the solver could become class (not instance) attributes.

But how do you run a model with a method of ModelParameters?
The same way we are already doing it from a function:
load the solver class and call complete on it.
But pass it the 'as_action' argument as well.

Wait! Uptodate needs to specifically exclude some of the slots.
We don't want the whole dictionary. Just parts of it.
So we'll need to specify which keys go into that one.
So, it does need to be a method of some kind after all.
For example, there could be a "read_keys" attribute, or something like that.
It could be a class attribute, I guess.
ParameterSet can check this to make sure the required keys
were defined by the YAML file,
and it can have a method to return a dictionary with only those keys.
For the required keys check, the class attribute is "\_required_attr".
For the config, the class attribute is "\_config_attr"

How does taskname work?
It can be a read-only property.
The default getter method will read the class attribute taskname_source,

I'm not really happy that useful must now import from doit.
It makes doit less optional than it was before.
I mitigated this by putting it in a try block.

How do inputfiles and outputfiles work?
Probably read-only properties again.
But no default getter method; it's entirely class-dependent.
But you have the input file (maybe) at instantiation.
So maybe it's not a read-only property.
It's an (instance) attribute.
It gets a slot.
Initialized at instantiation.

The inputfiles entry for the yaml input file is: 'initializer'

Note that if you give the base class mutable class attributes,
they are shared by all derived classes.
Unless, apparently, you use a metaclass.
But these shouldn't be class attributes anyway: they should be instance attributes.

How about this:
what if inputfiles and outputfiles weren't dictionaries,
but rather lists of attributes containing input and output file paths,
respectively.
These could be class attributes.
Hmm. Sometimes they might need to be  nested, rather than first-level attributes.
Except that the file it was loaded from isn't an attribute.
But it could be.

So now inputfiles and outputfiles become read-only properties.

The attribute for the parent file is: 'sourcefile'

Ah, no, still not right.
The properties read in from the yaml file are not full file paths.
We need a transformation function to locate those.
So it's going to be more complicated than just a list of attributes.
Unless, along with that, you had a list of folders.
But that's not the way the sourcefile comes in.
You could split it up, of course.
But then we have things like the otherfiles dictionary.
Converting that into the proper list of input files isn't trivial.
And, to avoid duplication, it must be done in a way that
the rest of the function can use the resulting information.
This is where the idea of a dictionary came from.

Maybe after any initialization,
we need a function that processes only required attributes
(and maybe not even all of those)
and populates the inputfiles and outputfiles accordingly.
Maybe it should just be a list?
No, then you get redundancy.
So here is the fundamental problem:
the file-related attributes we read in are relative paths.

Aha! You could have a class-attribute dictionary of the folders
those are relative to, by attribute name.
That works with the setup of having a list of attributes.

But now, what about paramgen's otherfiles?
For that, I added the \_more_*files attributes, which are direct lists of files.

I decided the input and output files should come from a simple list of attributes.
If the children have them, you should just ask the children to provide them,
by calling the appropriate method on them.
That means there is now a need for class attribute listing attributes containing children to be queried.


# Code/Misc

_TODO_ code from notebooks into modules
- time-domain Fickian?
- steady-state PNP?
- time-domain PNP w/o reactions?
- time-domain PNP with reactions

_FEATURE_ 2D mesh generation

Also note that body-centered.yaml says to use body-centered.geo.jinja2,
but really that already only works for body-cen2.yaml, because it has a Z5.

Eventually, we will need a test problem for this.
(Should probably use Fickian solver.)

And then we'll need a new gmsh .geo template file for parameter generation.
Currently, mesh.yaml.jinja2 is specific to the nanopore geometry.

_ISSUE_ have initial potential consistent with other initial conditions, including boundary conditions.
Tried solving Poisson by itself first, but couldn't get results into the mixed function space.

_ISSUE_ figure out a way to get t=0 into the same VTK file as the other timesteps

_FEATURE_ scaling factor on x and y values to do unit conversions

_FEATURE_ add metadata text to model plots
a text box consisting of strings using a template,
which is rendered using the info dictionary itself.
(So the template should use the names.)

Just like with a legend, the challenge will be where (and how) to locate this.

_TODO_ replace 'from folderstructure import \*'
This is challenging, because you have to search multiple files for the strings therein.

_TODO_ find a way to get coordinates of the surface normal used in a flux calculation
The notebook dated 2017-11-06 is where I was working on this before.

_TODO_ complete the "test" analysis, even without middle surface
That is, run a face-centered model somewhere, even without the middle surface.

_TODO_ wiki page on flux integration over internal boundary
Should probably document how to do external boundary as well, for comparison.

_TODO_ syncing/mirroring
- set up a scratch folder on /u1 (could use unison to sync it): will speed up runs on CP233
- holly and/or dlx (has to be rsync, unless we can build unison on the server)

For holly/dlx:
sync entire code folder?
That would mean not having copies of login scripts.
The goals are:
- everything on the remote machine has a local storage place
- avoid duplication: don't keep multiple copies of anything on the same computer
(it's ok to have scratch and home versions of something, but try to keep scratch clean)
- keep the same folder structure everywhere, as much as possible

This suggests the need for:
- a way of quickly copying to scratch, syncing, and removing scratch for each computer
- an established "one folder structure" for everything not machine-specific
  Or not. Maybe it is multiple structures, the roots of which may be at different locations on different machines.
  Kind of like what I was hoping to do with my flash drive before.
- a storage structure for machine-specific files (login scripts, etc.), and a way to update them into their machine-specific locations as well.

The real question is getting unison on the remote machines.
In fact, it needs to be the same version, even,
on any machine you want to sync with one of those.

After that, you start here:
https://www.cis.upenn.edu/~bcpierce/unison/download/releases/stable/unison-manual.html#remote

_FEATURE_ running on holly/dlx
It would also be nice if we had a way, on these systems,
to break the doit tasks into chunks that could be assigned to nodes.
You could do it by giving each node a different symlink for control, but
- you still need to modify dodo so that it can take an argument telling it where to find control.yaml (or other filename)
- it would be nice to have more granularity than that
- the database for each one has to be different, and then they have to be re-synced

See also `doit -n` to allow tasks to run in parallel.

dodo files can indeed take command-line parameters
http://pydoit.org/task_args.html

doit also has command line option `--db-file` to specify the database file.
There is a another option to specify the database format.
The default uses the python `dbm` module,
so if nothing else you could write your own code to reconstitute the new database.

You could certainly split one of the massive input yaml files at the time
of parameter generation.
Just do a round-robin on the output files.
The assumption here is that adjacent runs will have similar resource requirements.
That isn't always true, but it should give better load-balancing than splitting
the file after generation.

Note, however, that holly doesn't have many nodes, but they each show 64 cores.
And, of course, you can allegedly try mpirun for fenics (see item below).

That would require a change in how the fenics tasks are started.

_TODO_ use mpirun for fenics calculations?
(not for pre-and post-processing)

https://fenicsproject.org/qa/6969/how-to-go-parallel
https://fenicsproject.org/qa/8459/what-is-the-simplest-way-to-use-mpi
https://fenicsproject.org/qa/3025/solving-a-problem-in-parallel-with-mpi-&-petsc

_EFFORT_ we need to run for face-centered geometry as well
This requires adding the interior surface to this mesh,
which means redoing its geometry.
Is there a better way?

While I'm at it, I should redo the body-centered geometry to be more logical.
Maybe even include two additional internal surfaces,
at the ends of the pore(s).

Don't forget to remove mesh/geomdef/body-cen2.yaml and update
params/mesh/debug.yaml accordingly.

#Specific Equations
## Homogenized Fickian

_ACTION_ see discussion in Auriault
_ACTION_ look through homogmwe again as well

# Problem Description

_TODO_ Smol: use z for charge instead of q --wait! this conflicts with z as a coordinate! (ok in PNP b/c of species index)
_TODO_ reaction terms (not written up yet): find a good index name for reactions
_TODO_ would it be better to use delta_... for test functions (ala variational calculus) instead of v?

It has its own _TODO_ list.

## Physical Surface/Volume Figures

_TODO_ For the part where we have figures of the geometry,
it would be nice if these could be auto-generated from the lattice yaml file.
That way, if I add internal surfaces, the drawings could auto-update.
Of course, right now the yaml file doesn't contain the physical dimensions.

This would be a time-consuming effort.
Is there any other way to generate the necessary figures?
Obviously, you could do the same thing, but without automation.
The downside to that is that you have to redo it if something changes.

So the real question, then is can we use a less complicated figure somehow?
Or rather, how can we make the figure less complicated,
while still showing what it needs to show?

The ideal solution would be to use the gmsh model itself.
Is there a way to get it to generate desired views?
One issue is that it uses the sequential point numbers;
it doesn't know what I called them.
On top of that, getting it to show what you want seems tricky.
If you shade the surfaces, sometimes the surface numbers pass behind them,
or are truncated.
You can get all the surface the same color, by shading by type instead,
but that doesn't fix the problem with the numbers not always being visible.

So, the options are:
- Do all the needed geometry figures manually, and redo them if they change
- Do all the needed geometry figures automatically, which will take too long
- Do the surface figures automatically, but the 3D views manually.
- Find a decent 3D cad program that can display things nicely, which may not exist

Note that for automatic generation,
the macro used for results figures would solve the issue
of the figures potentially being missing.

## Results organization

Maybe the result figures and discussion should be a separate document.
That way, the problem description could be generated before the analyses are run.

You could even set it up so that the document that requires completed runs
also pulls in the problem description itself,
so in the end you do get a document that has everything.
This would require splitting the description up as well:
- the contents that are part of both documents
- the part that pulls that in as a standalone document

For now, I stored what little I had for that section in results.tex

Decision: Include them in the same document,
with a macro that uses placeholder figures for missing images.

There may be ways to improve the current macro.
For example, maybe there should be a message that explains
why the graphic is missing (need complete calculations)
References:
https://tex.stackexchange.com/questions/176201/how-to-create-placeholder-for-missing-figure
https://tex.stackexchange.com/questions/44195/placeholder-for-figure-includegraphics

Also, case-specific inputs should be presented in the relevant results section(s).



## Biblatex reference
https://www.sharelatex.com/learn/Bibliography_management_in_LaTeX
https://www.sharelatex.com/learn/Biblatex_bibliography_styles
https://www.sharelatex.com/learn/Biblatex_citation_styles

# Post-processing

_FEATURE_ More general function to extract data for 1D plot
- needs to return two 1D arrays (independent and dependent variables)
- could give it coordinate axis, limits, and number of points
- but more generally, could give it two points, and a number of points
- in some cases you may want to put this in a loop, to run multiple lines and plot them on the same set of axes
- could also want different components of a vector on the same set of axes
- need to restrict to points that are actually inside the mesh or on the boundary?

_FEATURE_
A similar thing would be nice for 2D slices.
But here, masking points outside the mesh would be even more important.
This can be done with `tree=mesh.bounding_box_tree()` and `tree.collides(Point(...))`.

More difficult, but also very helpful, would be finding the actual boundary points themselves.
How could we do that?
Here is a package that *might* work:
https://github.com/mikaem/fenicstools/wiki

Another idea would be to get the points from paraview.
Maybe it can get points on a cut-plane.
Or, maybe gmsh can do it. You only need the mesh, not the solution.
(See http://onelab.info/pipermail/gmsh/2008/003822.html:
it seems you can export plots in raw text.)
In any case, this is pretty low-priority right now.

# Mesh
_ACTION_ Look up "field" in gmsh tutorials to try to resolve issue with centerline. (tutorial 10, and the manual discussion on controlling mesh size)
_ACTION_ check for compatibility of gmsh version (eg Ruled Surface vs Surface)
_ACTION_ add validation of geometric inputs

# Someday/maybe

_maybe_ solve with flux as another primary variable
_maybe_ calculate potential for Smoluchowski using nonlinear Poisson-Boltzmann

_maybe_ doctests? some other kind of test?
Maybe just create scripts to do the basic validation on buildgeom and the solvers.

_maybe_ mesh refinement study
_maybe_ study of required H value

_maybe_ Find a replacement for codenamize
I experimented in the notebook with some other stuff.

_maybe_ investigate more accurate time-step methods, such as Adams-Moulton

_maybe_ look into opencascade and gmsh for mesh generation
https://en.wikipedia.org/wiki/Open_Cascade_Technology

opencascade is not a program, but rather an SDK.
gmsh refers to it as a geometry kernel.

Using this kernel, gmsh can import BREP, IGES and STEP files (I think).

There is a set of python bindings called Python OCC:
http://www.pythonocc.org/
Precompiled binaries are available from conda-forge.

There is a "scripting interface" called pycado,
but it looks unfinished
http://julienbld.github.io/pycado/

FreeCAD is a program that can generate IGES and STEP files.
It is scriptable with python,
and uses the opencascade kernel.

It might even be possible to get something out of Blender,
but maybe not.
https://docs.blender.org/manual/en/dev/data_system/files/import_export.html

--------------------------------------------------------------------------------
# Non-active items only below this point

_question_ is it possible for fenics to give us the condition number of a system?
If nothing else, can we get the relevant matrix and calculate it from that?
