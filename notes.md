
__TODO__ scripts or models to generate the spatial_D output file
also, move it to a more appropriate location

_TODO_ add number of dimensions to geometry defnition file, and from there to mesh metadata

_TODO_ switch to ruamel.yaml, and update the wiki

_TOOD_ zipping large data sets may need to come under control of doit
Or should it be available in the simulator?
What if we're not using doit?
The trouble is it's a command-line thing, not python.

_TODO_ time-selection wrapper for datasteps (see below)

_TODO_ change species_info to a list of species dictionaries, which become species objects

_TODO_ see note in output_eff.fluxfield

_TODO_ command line "select" argument: apply more directly - don't instantiate objects first
This is actually not easy to do.

_ISSUE_ the data folder structure, and how different attributes specify different parts of it, can be confusing
- The input yaml file's own name defines the basename.
  Maybe that feature should be removed, requiring the basename to be specified.
  And call it something more specific, like "basefolder".
  Actually, the basename could be included in the yaml file already, but:
    - currently, it would be overwritten by the file's basename
    - that doesn't help with the case where the model and mesh have different basenames
        This suggests just providing the mesh name and directory.
- "meshname" is used to compute filenames.
  Maybe instead the file names themselves should be listed.
  I already had some discussion about this under the notes on pathlib.Path below.
The fundamental question is:
Should we get rid of folderstructure, and require the input files to use hard-coded paths instead?
Yes, this will increase the number of input parameters, by making everything more explicit.
It could also make the input files more machine-dependent.
Unless we keep the datafolder, which can come from an environment variable if necessary,
and make all paths relative to that.

_TODO_ rename paramgen to genparams everywhere
Not just directory and filenames, but in the code itself.

# Formula derivations
- _TODO_ NP linearization notebook
- _TODO_ start reaction rate function linearization notebook


# Equation Builder & Loading States

Rather than having separate codes with hard-wired weak forms,
maybe we want to be able to select weak form components.
The PNP linearization notebook (which otherwise doesn't work yet)
does some of this.
To solve Smoluchowski instead of PNP, just solve an initial potential,
then no longer include the Poisson equation.
(And more generally, establishing initial conditions should always happen for time-domain problems,
  and for iterated attempts at steady-state.)
For Fick's law, omit the terms with charge.
(The TDPNP simulation now does this for cases of zero charge.
  We could do it better, though.)

Or maybe we allow the 'equation' argument to be more than just a string.
Maybe that's where we specify what is and is not included in the equation.
Maybe there's even a way to allow customizations to contribute equation terms.

Or maybe one allowed "equation" is 'builder',
which allows the equation terms to be specified in the conditions field.
Maybe "equation" should be changed to 'sim_module' for clarity.

In fact, more generally, we'd like to be able to build models up from commands,
similar to the way the data extraction takes place.
That makes the input files something like a DSL.
Which ultimately means what we want is a library of tools,
which we can call from a python script.
So we're just creating abstractions on top of FEniCS.
And yet, we do still need data files.
Where do we draw the line between code and data?

To some extent, that's what the simulator modules already do.
And it's what the most recent notebooks do as well.
For right now, just focus on tools to make doing those things easier.

This will also require doing some verification of the input.
Check to make sure there are sufficient equations for each unknown, etc.
How could that be done?
For example, anything that doesn't have a diffusion constant must be involved in a reaction equation.
And vice-versa.

Maybe as a first step, we just need to make it easier for simulator modules to construct weak forms.
Toward this end, a "library" of weak forms might be helpful.
But they'd need to be functions, instead of UFL forms,
to be passed in arguments such as the TrialFunction or Function, and TestFunction.

The goal is to reduce redundancy.
If I have a steady-state Smoluchowski equation solver, and a time-domain one,
and then the same for PNP, all of them need the Smoluchowski weak forms.
So maybe we want functions that return the terms for an entire equation.
No, then time-domain and steady-state would be different.
We want what I've called the steady-state weak forms.
Steady-state simulators just set this equal to zero in their equation.

Once we have such a library (module),
the next thing we need is an attribute of ModelParameters that specifies which forms to include.
We also need a good way to support construction of zeros for the linear forms.


So these two things need to work together: the weak forms library, and the way settings are listed in ModelParameters.
At top level, ModelParameters will have a new attribute "weakforms".
This will be a sequence, with each item specifying a term, or group of terms,
to add to the weak form.
At the lowest level, the library will contain functions to return terms or groups of terms,
with arguments providing the actual TrialFunction or Function instances and TestFunction instances.
In between, we need:
- A way to construct entire equations, so that we know what zeros are needed in the linear form
- A way to specify attributes to be used in the function arguments
Maybe the list in ModelParameters really is an equation list,
and each equation in turn calls the necessary functions for its individual terms.
That way, reusable equations can be built up from reusable parts.

# Fick's Law with Reactions
- DONE: set up an example problem (in debug, I guess, base it on debug03)
- DONE: only 2 species: Ca and CaCaM
- DONE: reaction rate code based on the one for ill-sleep
- DONE: copy the latest PNP linearization notebook
- work out correct weak form
- revise the weak form to take out all PNP
- I think I did a TD Fick's Law simulation in a notebook once, so you could borrow from that
- get the reaction stuff from the TDPNP module
- switch to nonlinear problem and solver
- once it is working, create a module for it
- then finish the documentation

Or do we want the equation builder instead?

# validation tests
set up problems on the unit square with known solutions:
- fick's law: such that correct solution is a simple gradient
- smoluchowski: such that correct solution is a Maxwell-Boltzmann distribution
- PNP: probably need to set up a problem that can be solved analytically, and do so. Or is there one already published?
- hom. fick: use a single mesh from exotic-earth, or maybe the whole thing?

What about reactions?
Have a reaction-only simulation, without diffusion?
It depends on what your simulators are capable of.
Maybe a Fickian reaction-diffusion simulation could be contrived to be solvable.
Beyond that, maybe not.

Don't overwrite debug, but rather create a new stem: validation

_TODO_ Need a script to compare output to expected values
`validate.py`
Rather than storing the validation data in the script,
store it in a dedicated (tracked) location,
and have the script just make the comparisons.
This may require some post-processing of the model runs,
to generate a comparable data file.
And potentially, then, some separate scripts to generate the data files with the expected results.
So, maybe we need a separate "validation" directory in src
to contain these scripts (and the comparison script),
and a "validation" directory in data to hold the expected results,
for comparison to what will be generated in data/postproc/validation.
Not in src: we probably just need a "scripts" directory in general,
for other cases where I have to generate data.

# anisotropy
- fickian_unhomog
- smol_unhomog
- fickian_homog
On hold for now per Pete's request.

Accept diffusion constants in three forms:
- single value for isotropic conditions
- list of values (length matching dimensions) for diagonal matrix (order x,y,z)
- full matrix (user must do any coordinate transformations themselves)

# homogenization (exotic-earth)

_TODO_ do data extraction so we can actually use the results
_TODO_ then do data collection and plotting

Then, we need to try to abstract/generalize this more, somehow.
- can we set up a general structure for rectangular periodic boundaries?
- spatial variation in diffusion matrices? or properties that vary by cell as a generalization of that?

_TODO_ set up as a validation test (see above)

# ill-sleep (and debug of the PNP-reaction simulator)

_TODO_ compute integrated total charge for each step
For problems with reactions or neumann BCs on the species this might not be conserved,
but it could still provide insight.

_TODO_ linearization
What is the argument structure now for the linearized reaction rate functions?
They need the whole vector at both steps.

_TODO_ allow BT to be an Expression instead of just a float.
More generally, this is a candidate for refactoring:
write a function that does the check for parameters as lists and returns an Expression or a float accordingly.
It will need access to the element. Or can it take a FunctionSpace instead?
Either way, that means it will probably need to be a method of the simulator.
Actually, let's just pass the element, as there could be times we need different elements in different cases.
It still needs to be a simulator method, though, because
it will require a more general way to track the expression arguments,
and make sure they are updated as needed.

Argument tracking:
all_expressions = {expression: arguments}
expression_argument_values = {argument: value}

_TODO_ refactor out re-usable components
Weak form portions that could be used by various PNP simulators:
Imagine we have multiple ways of solving PNP, and they have a lot of basics in common.
This may need to wait until that's actually the case, as that will make it more clear what's reusable and what isn't.

Maybe some of this is the equation builder instead.

Neumann and Dirichlet boundary conditions:
  - store in attribute, or return value?
  - take argument, or get directly from conditions?
  - how to handle mixed and non-mixed elements?


_FEATURE_ exponential time steps
But dt is in the weak form.
So if it changes, does the weak form need to recompile in FFC?
If so, that would argue for stepwise variations.
But what if instead we put explicit times into the weak form,
rather than just dt?
Actually, it seems to work with dt as an Expression.
Well, it worked in a simple example I did in a notebook.
When I tried to put it into PNP, it didn't converge,
even for a single step.
That branch of the code still exists,
so you can look back at it if necessary.

Pete did this in `/net/share/shared/labscripts/dolfin/template_adaptivetime.py`
using a `Constant` as opposed to the `Expression` I attempted it with.
Maybe that makes a difference?

New idea: just put it in the weak form explicitly.
That is, transform from t to tau by an exponential function,
and get a weak form that only has tau.
Then take constant steps in tau.
(See notes in notebook dated 12-Mar-2018.)

_TODO_ for debug, consider using an expression to specify the initial electric potential.
Requires code change to allow Expressions for Dirichlet conditions just like they are for Neumann.

_TODO_ automatic color selection for the plots with lots of time-steps as different series

_TODO_ time-selection wrapper for datasteps
The catch is how to specify those times.
(Because doing all of them is too much.)
You need to specify some attribute that will be used.
You could try just a list of values for that attribute.
But in some cases that might mean a very long list.
And, if the values are floats, they'll never match exactly.
So that approach could work for short sets of integer (or string) values.
But for time steps, we need something else.

You could keep track of a "next snapshot time",
and whenever the current value is found greater than that one,
you take a snapshot, and then update the next time.
So the specification needs to include the delta.
Or, a list of floats could work this way, as long as it wasn't too long.

Basically, this issue of selecting steps for output, rather than doing every one of them,
could be a wrapper around any of the time-dependent stuff.
So maybe what we need is a wrapper function,
which defines when to generate output,
and then accepts a command list.
Or maybe datasteps itself is that wrapper.
It takes a series of pairs.
Item 0 is the definition of when to perform the list of commands in item 1.

Maybe this could work something like stopping criteria.
You provide trigger values, and step values.
The triggers a just like stopping criteria.
But after each trigger, the step is added to get the next trigger.

See the `stopnow` simulator method and `StoppingCriterion` class in tdpnp_unhomog,
and also the enhancements in a notebook such as `p20180607_SS_PNP_opsplit`.

_TODO_ This needs to be documented somewhere:
Diffusion exclusions:
just put "null" in as the diffusion constant.
pyyaml will convert that to a None.
Then just exclude such terms from your weak form.

_TODO_ Do a check to make sure that the simulation is electrically neutral at first:
take sum of initconc*z.
Maybe you do this in your input generation?

_FEATURE_ Time-dependent Neumann condition on select surfaces:
Expressions: https://fenicsproject.org/pub/tutorial/html/._ftut1006.html#ch:fundamentals:diffusion
simpler expressions: https://fenicsproject.org/pub/tutorial/html/._ftut1011.html#ch:poisson0:DN
restrict surface: https://fenicsproject.org/pub/tutorial/html/._ftut1014.html#ch:poisson0:multi:bc

It would be easy to just check for a string, and convert it to an expression.
The problem is that if it is an Expression, you need to give it parameters (t in this case).
And those parameters then have to be updated in the expression object at each step.
Maybe that could be done with a custom function in datasteps?
We could make the condition itself a tuple instead of just a string:
- the first argument would be the expression string
- the second could be a dictionary, sent as kwargs

# Code/Misc

_TODO_ more general reaction functions approach
- simulations with different species overall may still contain the same reactions
- sometimes we need to distinguish between previous and current step, but not always.
Ideas:
- use a dictionary? Seems excessive, and doesn't resolve the issue of multiple steps.
- specify a mapping to the arguments. But what is that in terms of?
Multiple steps seems like it will always require a separate function set up for such.
So maybe we just want a way to keep straight which species is which.
We could have the mapping to argument locations specified in the input file.

_TODO_ Separate "model" parameters and "simulation" parameters
For example, a simulation timestep would be a simulation parameter.
A boundary condition would be a model parameter.
Maybe time-domain or steady-state is even a simulation parameter rather than a model parameter.
(Although that could get tricky, with time-dependent boundary conditions.)

_TODO_ The ability to resume time-step jobs would be nice,
but this requires being able to save out the complete state.
Or maybe just be able to load the uknowns, I guess.
You can reload the mesh and BCs, etc, from scratch.
In fact, based on Pete's suggestions, we might even switch equations.
But somehow, you need enough values to get back the original functions.
Would it work to just store the nodal values?
I'd have to do a test to see.

_TODO_ subclass of pathlib.Path for paths

p=pathlib.Path(...)
- folder: p.parent
    Except that if p is itself a folder rather than a file, this will give the parent folder, not the same folder
    You could check to see if it is a folder, but that requires checking the disk.
    Your subclass could have a boolean flag for whether it is a file or folder,
    which you'd have to maintain at each operation.
    Or, you could just keep track of it yourself.
    That is, know what is supposed to be a folder and don't call folder on it.
    Maybe just checking the disk is the best approach.
    But what if the path doesn't actually exist yet?
    Then checking the disk won't tell you anything.
    (See below for continued discussion.)
- stem name: p.stem (note that for multiple extensions this only takes off the last one)
    If subclassing, you could create a read-only property that uses p.stem.split('.')[0] to get the base name
- extension: p.suffix (note that for multiple extensions this only gives the last one, and it does include the dot)
    If subclassing, you could create a read-only property that calculates this
- file name: p.name0
    Of course, for a folder, this isn't a file.
- full path: str(p) (repr(p) is different, and should not be used)
- assure_dir: create the necessary directories if they don't already exist.
    And this means having disk operations is ok.
    It also requires knowing if the endpoint is a file or folder.
    (See below.)

You have to use this class as the appropriate object attribute,
meaning at initialization it has to be there.
That doesn't fit with the way ParameterSet works.
Although we could read in the string, then at initialization change it to the class.
Or create a different attribute for the class.

Towards that end (initialization), maybe the input files should be more consistent.
Sometimes they give stem names only, and sometimes they give filenames.
For paramgen, you even need to specify the folder (since it doesn't know what type of parameter set you want to make otherwise).
A good compromise would be that the extension should always be included.
That is, it's always a path relative to the respective location given by folderstructure.
In many cases, that means all it is is a filename.
But in the case of paramgen, it needs more than that.
Alternatively, you could add a folder field to paramgen input files, but that seems unnecessary.
Another challenge is that there are things like meshnames which do get extensions added sometimes, but not always. 
So maybe we just do what we're already doing, but with Path instead of with just strings.

Then, at object initialization, those filename parameters are converted to Paths,
by appending the folder specified in `_folders`.
But sometimes you need to add an extension as well.
(This could even become an optional step in ParameterSet init: skip it for empty or nonexistent entries.)
Instead of calling a function to get the full path, you use the appropriate attribute of the Path object
or its subclass.
Rather than converting, what if we add a `_files` or `_fpaths` attribute.
The original entries can remain, and this attribute is a dict (or namespace)
with fpath versions of the files.
That's great for direct attributes.
But what about `_more_inputfiles` and the like?
Those are required to be complete paths (as strings),
so they'll be refactored on a case-by-case basis.
So the new method can basically work like common.ParameterSet.full_path does now.
That is, for a given attribute list (which is usually stored in an attribute),
it prepends the matching entry in `_folders`, if any.

Everything in folderstructure should become a Path (or the subclass).
If it's going to be the subclass, then folderstructure would have to import common.
Which means it needs to locate the source directory first.

Implementation steps:
- DONE: new depcheck on pathlib
- DONE: new class based on pathlib
- use in folderstructure
- method of ParameterSet to initialize objects from their input strings
- call that method in all subclasses of ParameterSet
- revise docstring of ParameterSet.file_list
- replace calls to ParameterSet full_path with appropriate object method
- delete ParameterSet.full_path
- for all classes derived from ParameterSet, look at all input/output file attributes, and catch each use, to use the correct object method
- check all uses of os and osp (definitely os.makedirs should go)
- remove unneeded imports of os and osp

_TODO_ should common be split up into more than one module, possibly inside a package called `common`?

_TODO_ should ParameterSet be split into a base class,
and a derived class that includes all the doit support?

_TODO_ the LPB simulator dataextraction output files don't get listed as targets currently.
Search conditions (actually, it's immediate children) for dataextraction.

_TODO_ debug01 postproc needs to make PDF as well as PNG

_TODO_ code from notebooks into modules
- time-domain Fickian?
- steady-state PNP?
- time-domain PNP w/o reactions?

_FEATURE_ 2D mesh generation

Also note that body-centered.yaml says to use body-centered.geo.jinja2,
but really that already only works for body-cen2.yaml, because it has a Z5.

We'll need a new gmsh .geo template file for parameter generation.
Currently, mesh.yaml.jinja2 is specific to the nanopore geometry.

_ISSUE_ have initial potential consistent with other initial conditions, including boundary conditions.
Tried solving Poisson by itself first, but couldn't get results into the mixed function space.

_FEATURE_ scaling factor on x and y values to do unit conversions

_FEATURE_ add metadata text to model plots
a text box consisting of strings using a template,
which is rendered using the info dictionary itself.
(So the template should use the names.)

Just like with a legend, the challenge will be where (and how) to locate this.

_TODO_ find a way to get coordinates of the surface normal used in a flux calculation
The notebook dated 2017-11-06 is where I was working on this before.

_TODO_ complete the "test" analysis, even without middle surface
That is, run a face-centered model somewhere, even without the middle surface.
Or, get rid of it.

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

_FEATURE_ use mpirun for fenics calculations
(not for pre-and post-processing)

https://fenicsproject.org/qa/6969/how-to-go-parallel
https://fenicsproject.org/qa/8459/what-is-the-simplest-way-to-use-mpi
https://fenicsproject.org/qa/3025/solving-a-problem-in-parallel-with-mpi-&-petsc
https://www.allanswered.com/post/eqnvq/get-a-displacement-at-a-node-in-parallel/

https://fenicsproject.org/olddocs/dolfin/1.6.0/python/programmers-reference/cpp/common/MPI.html

The challenge is that the python script will be run in parallel,
meaning some operations that are only intended to happen once,
happen multiple times instead.

Pete mentioned on Slack that he has an example code that does this correctly.
(Thursday, February 15, between 12pm and 12:30 pm)
He also said it's in `template_timedep.py`.

Got things working to generate meshes in HDF5 format,
and to have the simulator read those in.
But now the issue is data extraction: each process only has part of the mesh.
See log 2018-05-29.md

So, what will it take to get this working:
- install fenicstools
- use fenicstools Probes to extract data
- MPI gather results from different processes (or does Probes do this for you?)
- only rank 0 process should write output file

Alteneratively, I could try to see how Probes in fenicstools works,
and write similar python code.

_EFFORT_ we need to run for face-centered geometry as well
This requires adding the interior surface to this mesh,
which means redoing its geometry.
Is there a better way?

While I'm at it, I should redo the body-centered geometry to be more logical.
Maybe even include two additional internal surfaces,
at the ends of the pore(s).

Don't forget to remove mesh/geomdef/body-cen2.yaml and update
params/mesh/debug.yaml accordingly.

_FEATURE_ a new class like ParameterSet,
but not using slots.
Instead, use a jsonschema to indicate what keys are required on load,
and what structure is allowed.
http://json-schema.org/
https://github.com/Julian/jsonschema
(already present in anaconda)

The schema for each class would be a class attribute.
In init, run the validation on the input dictionary before assigning its values to the attributes.

Use mro to assemble schema from all parent class schemas (and the class itself).
The complete schema would thus be a read-only property.
Or is there a way to call the validation methods of all supers?
Probably not, as you'd need to know what keys are allowed by other classes in the mro.

Or is it better to stick with what we have now?
Slots might actually be more portable than using jsonschema.

So what would be the advantage, then?
The validation would be more thorough than just listing required keys and allowed keys.
It would also include type information. (And maybe more.)
Also, because there are no slots, you could still add attributes at run time.
The validation would only happen at instantiation, to confirm the data read in from the file.
(Towards that end, maybe you want a validation method you can run later on as well?
It would require converting the object to dictionary and validating that.)
I'm not sure that's enough of an advantage.

Maybe the class docstring could be yaml (with comments for actual text)
that includes the schema, to reduce redundancy.
Now that might be enough of an advantage.

Keep in mind that you won't be able to do 2nd-level inheritance very well this way.
That is, your base class will be a kind of meta-class.
(Maybe it should be exactly a meta class?)
But once you have a class built from a schema,
there's no point having another class inherit from it.
It can't really add to the schema.
(Unless, maybe the metaclass teaches it how to?)

So, a metaclass?
The `__new__` method of the metaclass will accept a string
which becomes the docstring of the class,
and is also used to construct its schema.
Each instance of the metaclass (the schema-based classes)
have a method that can check their current state against the schema.
This method should be called at the end of their init.
You would want them to inherit from ParameterSet or something like it, still.

_MAYBE_ can we come up with a better name for paramgen_tmpl?
- 'templates' is too generic: these are specifically for paramgen
- but we don't want anything that starts with 'paramgen' as its too close to 'params'

_MAYBE_ maybe EquationTerm, rather than having boolean attributes, should just have an attribute sequence.
Individual attributes can be added to this sequence.
You may want to set up a singleton that defines the attributes you want, so you don't always have to quote them.
Selection then just looks for those attributes in the sequence.
I suppose it does make nonsensical states more possible:
eg, you could accidentally set both 'linear' and 'bilinear' attributes on a term.
So, what's the advantage over the current approach:
no kwargs, so no slots, no subclassing
But then we do need the attribute class.
Also, right now we can do more than just booleans with the attributes.
To some extent, they can do what I was thinking about doing with the names.
Conclusion: leave it as-is for now.

# Problem Description

We need to reconsider this whole approach to the document.
The various goals are:
- document the math behind the code at a high level
- document extensive formula derivations (in a appendix, referenced in the Part above)
- document the code itself
- document the specific analyses performed with the code, including possibly validation problems

We really want this to be modular.
Right now it's kind of all mixed together.

Documentation of specific analyses:
- geometry
- mesh (see below)
- specifics of general equations from the math derivation (e.g. specific reaction rate functions)
- input values of constants, in some cases (e.g. eps_r, beta or T, etc.)
- expected results
- actual results (which require that the model has been run)

So how do we modularize weak forms, time-domain calculations, etc.
For example, sometimes we have special cases:
- steady-state instead of time-domain
- isotropic instead of general diffusion matrix
- linearization of a nonlinear equation

So maybe we start with a general equation, and then do special cases?
Really, it's more like we want to specify re-usable elements of equations.
For example, I only want to derive the PNP weak form from the PDE once.
I want to be able to re-use that same weak form in steady-state problems,
and time-domain problems.
I may do different linearizations on it.
So where does a particular linearization of a particular weak form,
particularly steady-state or with a particular time discretization go?
We need general forms of all these things,
which those special cases then reference.
So in terms of order, we need the most general concepts first,
followed by special cases.
That is, we need general discussion of linearization
before we linearize a particular equation.

_TODO_ bring in references from Zotero.

It has its own _TODO_ list.

_MAYBE_ Can we use a two-column format for the derivation portions?
Left side (narrower) is text, with equations on the right.
Maybe an invisible table?
How would this work with equation numbering?
Each row of the table should be numbered.
So could we get the equation numbers in their own column?

_MAYBE_ Smol: use z for charge instead of q --wait! this conflicts with z as a coordinate! (ok in PNP b/c of species index)
_MAYBE_ would it be better to use delta_... for test functions (ala variational calculus) instead of v?

Sphinx tips:
- docstrings in markdown format? the other documents? toctree?
  http://www.sphinx-doc.org/en/master/usage/markdown.html
  https://github.com/jxltom/sphinx-markdown-extension
- include latex files (only for latex output): use `raw`
  https://stackoverflow.com/questions/45064321/sphinx-including-tex-file-via-raw-latex?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
- cross-references
  https://stackoverflow.com/questions/21289806/link-to-class-method-in-python-docstring?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
- Changing the top level of the toctree to "parts" in latex.
  https://github.com/sphinx-doc/sphinx/pull/4733
- Citing references
  https://sphinxcontrib-bibtex.readthedocs.io
- genindex and modindex
  (only slightly helpful)
  https://stackoverflow.com/questions/36235578/how-can-i-include-the-genindex-in-a-sphinx-toc?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
  https://stackoverflow.com/questions/25243482/how-to-add-sphinx-generated-index-to-the-sidebar-when-using-read-the-docs-theme

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

Use tikz for generating 2D figures.
There are plenty of examples using it for 3D drawings as well.

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

So, those aren't really "results" subsections,
as much as they are specifics of solved problems.

I really think this is needed now,
but I wonder if these should be promoted to Sections.
Especially since now we have Time Domain and Reaction-Diffusion problems.
They can reference all the preceding Sections they need to,
including the relevant geometry.

So the way to think about this now is that we have general
sections that could apply to multiple problems,
and therefore do not contain results.
But then we have other sections where specific problems are described,
by references to the relevant sections.
Some additional inputs are listed,
and in some cases the specific weak-forms need to be shown as well.

_TODO_ set up boxes in the problem description, like I did in the Iridates calculation.

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
(Looks like the new gmsh version resolved this.) Look up "field" in gmsh tutorials to try to resolve issue with centerline. (tutorial 10, and the manual discussion on controlling mesh size)
_ACTION_ check for compatibility of gmsh version (eg Ruled Surface vs Surface)

# Someday/maybe

_maybe_ solve with flux as another primary variable
_maybe_ calculate potential for Smoluchowski using nonlinear Poisson-Boltzmann

_maybe_ doctests? some other kind of test?
Maybe just create scripts to do the basic validation on buildgeom and the simulators.

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
