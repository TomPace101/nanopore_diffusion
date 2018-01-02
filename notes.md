
# Code/Misc

_TODO_ code from notebooks into modules
- time-domain Fickian?
- steady-state PNP?
- time-domain PNP w/o reactions?
- time-domain PNP with reactions

_FEATURE_ 2D mesh generation

_ISSUE_ have initial potential consistent with other initial conditions, including boundary conditions.
Tried solving Poisson by itself first, but couldn't get results into the mixed function space.

_ISSUE_ figure out a way to get t=0 into the same VTK file as the other timesteps

_FEATURE_ investigate more accurate time-step methods, such as Adams-Moulton

_FEATURE_ scaling factor on x and y values to do unit conversions

_FEATURE_ add metadata text to model plots
a text box consisting of strings using a template,
which is rendered using the info dictionary itself.
(So the template should use the names.)

Just like with a legend, the challenge will be where (and how) to locate this.

_ISSUE_ doit tasks can't generate meshes when there are no models
This was kind of as intended: it only generates meshes that are needed.
But maybe you should generate all the meshes for the basenames requested in control.yaml.

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

_TODO_ find a way to get coordinates of the surface normal used in a flux calculation
The notebook dated 2017-11-06 is where I was working on this before.

_TODO_ complete the "test" analysis, even without middle surface

_TODO_ shell command files for generating .msh and .xml files, which doit then calls?

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

_TODO_ look into opencascade and gmsh for mesh generation

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

_TODO_ Smol: use z for charge instead of q

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

_mabye_ Find a replacement for codenamize
I experimented in the notebook with some other stuff.


--------------------------------------------------------------------------------
# Non-active items only below this point
