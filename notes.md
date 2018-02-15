
_TODO_ switch to ruamel.yaml, and update the wiki

_TODO_ set up boxes in the problem description, like I did in the Iridates calculation.

# ill-sleep (and debug of the PNP-reaction solver)
Other stuff:
- I need to figure out how to exclude CaCaM in the weak form for diffusion but not electric potential
- then I need to figure out how to specify the expression for the Neumann boundary condition.
- and at some point I need to fix the weak form based on my earlier observations
- and I need to set up the necessary post-processing routines, to generate plots (only model plots in this case, I think)

Diffusion exclusions:
just put "null" in as the diffusion constant.
pyyaml will convert that to a None.
Then just exclude such terms from your weak form.

Do a check to make sure that the simulation is electrically neutral at first:
take sum of initconc*z.
Maybe you do this in your input generation?


# Code/Misc

_FEATURE_ parametric definition of mesh locations

There is some awkwardness in that buildgeom needs to put the location of this file into the .geo file.
So both buildgeom and geom_mk_mesh have to calculate the location.
This also makes the .geo file machine-dependent.
If you switch to a different computer or move the whole project,
the .geo files have to be regenerated.
Maybe gmsh could accept a relative path?
ie "../../paramlocs/{basename}/{meshname}.yaml"
__TODO__

The existing profile outputs should be refactored to use this. __TODO__

_TODO_ refactor solver modules into a package?
That would group together the various modules in a more logical way.
But what complications would it cause?
The package should probably be called "solvers".
rxn_rate_funcs would need to go in there as well.

_FEATURE_ greater modularity in data extraction functions

Right now, everything is going into solver_general.
But not all the output functions there are appropriate for that, really.
This might be a good use case for multiple inheritance.
Or, maybe even better, have a module of just extraction functions, which are not part of a class.
Then, in the class definitions, you can just assign them to class members.

The reason for this is that the extraction functions not only depend on the equation,
but also in some cases on the geometry of the problem as well.
From an inheritance perspective, there are base classes appropriate for an equation,
and then derived classes with data extraction methods appropriate to both the equation and the geometry definition (not the parameter values).
The use of parametric locations could help with this.

This could be handled by customizations.

_TODO_ use pathlib.Path for paths
Or maybe subclass it.

p=pathlib.Path(...)
- folder: p.parent
    Except that if p is itself a folder rather than a file, this will give the parent folder, not the same folder
    You could check to see if it is a folder, but that requires checking the disk.
    Your subclass could have a boolean flag for whether it is a file or folder,
    which you'd have to maintain at each operation.
- base name: p.stem (note that for multiple extensions this only takes off the last one)
    If subclassing, you could create a read-only property that uses p.stem.split('.')[0] to get the base name
- extension: p.suffix (note that for multiple extensions this only gives the last one, and it does include the dot)
    If subclassing, you could create a read-only property that calculates this
- file name: p.name
- full path: str(p) (repr(p) is different, and should not be used)
In fact, any function that needs a string instead of a path object will need
to get str(p.whatever) instead of just p.whatever.
Paths are immutable, but can be concatenated with "/"
and can get suffixes added wtih p.with_suffix()

You have to use this class as the appropriate object attribute,
meaning at initialization it has to be there.
That doesn't fit with the way ParameterSet works.
Although we could read in the string, then at initialization change it to the class.
Or create a different attribute for the class.

Conclusion: better than the way I'm doing it now, and part of the standard library.
This is what I should be doing.

Everything in folderstructure should become a Path (or the subclass).

Towards that end, maybe the input files should be more consistent.
Sometimes they give stem names only, and sometimes they give filenames.
For paramgen, you even need to specify the folder (since it doesn't know what type of parameter set you want to make otherwise).
A good compromise would be that the extension should always be included.
That is, it's always a path relative to the respective location given by folderstructure.
In many cases, that means all it is is a filename.
But in the case of paramgen, it needs more than that.
Alternatively, you could add a folder field to paramgen input files, but that seems unnecessary

Then, at object initialization, those filename parameters are converted to Paths,
by appending the folder specified in `_folders`.
(This could even become an optional step in ParameterSet init: skip it for empty or nonexistent entries.)
Instead of calling a function to get the full path, you use the appropriate attribute of the Path object.

The catch is that sometimes you really do want a string, and you have to remember to cast it.


_TODO_ should ParameterSet be split into a base class,
and a derived class that includes all the doit support?

_TODO_ in buildgeom, validate geometric inputs (different formulas for different geometries)

_TODO_ the LPB solver dataextraction output files don't get listed as targets currently.
Search conditions (actually, it's immediate children) for dataextraction.

_TODO_ the name of info.yaml appears in several places

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
He also said it's in `template_timedep.py`

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
