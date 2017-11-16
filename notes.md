
# Code/Misc

_ISSUE_ data/paramgen/model_smol_uh.yaml.jinja2 has some body-centered mesh-specific information in its bcdict.

_ISSUE_ doit tasks can't generate meshes when there are no models

_TODO_ units analysis as described below, to get an appropriate value of beta.
Then re-run, then generate the figure.

_TODO_ learn how &id001 works in yaml.
In particular, can it only refer to something in the same document?
If not, I'd like to use it to reduce redundancy in some of the yaml parameter files.
Particularly things like data extraction.

_ISSUE_ config_changed doesn't seem to be working for the plots.
Try it yourself and see.

_TODO_ add pore limits as vlines on all centerline plots
_TODO_ add pore radius as vline on all radial plots

_TODO_ plotting
Some plotfuncs will need additional info.
Put it in a dictionary?

maybe plotfuncs could be a sequence as well
call all of them before generating output

Yes, this should work more like the way data extraction does.
Plotting consists of a series of tasks.
In fact, maybe we should have a way to call the methods of the axes object,
and pass them a keyword argument dictionary.

_TODO_ plotdata.ModelPlotFigure.plot_radial_potential has some hardcoded numbers specific to the current body-centered mesh
A more general way to add hlines and vlines would help,
but not resolve the issue completely,
as those values need to come from a calculation.

Maybe model plots needs calcfuncs too, then?
How would that work, since you don't have a DataFrame?
As stated below, you'd like to be able to add new series as well,
not just hlines and vlines.

Maybe there should be hlines and vlines properties,
just like series.
And the calcfuncs can access them.

_TODO_ collection plot: what if I want to add a series?
This should probably be an option for PlotFigure in general,
as hlines and vlines

_FEATURE_ doit tasks for parameter generation
tasks are generated based on reading the output of such a task.

doit has a way to resolve this, of course:
http://pydoit.org/task_creation.html#delayed-task-creation

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

_TODO_ actually generate the plots we extracted data for (concentration vs centerline)

_TODO_ shell command files for generating .msh and .xml files, which doit then calls?

_TODO_ wiki page on flux integration over internal boundary
Should probably document how to do external boundary as well, for comparison.

_TODO_ syncing/mirroring
- holly and/or dlx (has to be rsync, unless we can build unison on the server)

_TODO_ look into opencascade and gmsh for mesh generation

_TODO_ use mpirun for fenics calculations?
(not for pre-and post-processing)

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
## Fixed potential

_ACTION_ get boltzmann constant into module?

## Homogenized Fickian

_ACTION_ see discussion in Auriault
_ACTION_ look through homogmwe again as well

# Problem Description

It has its own _TODO_ list.

_TODO_ work through details of consistent set of units, with lengths in nanometers.

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

_maybe_ maybe the result figures and discussion should be a separate document.
That way, the problem description could be generated before the analyses are run.

You could even set it up so that the document that requires completed runs
also pulls in the problem description itself,
so in the end you do get a document that has everything.
This would require splitting the description up as well:
- the contents that are part of both documents
- the part that pulls that in as a standalone document

For now, I stored what little I had for that section in results.tex


biblatex reference:
https://www.sharelatex.com/learn/Bibliography_management_in_LaTeX
https://www.sharelatex.com/learn/Biblatex_bibliography_styles
https://www.sharelatex.com/learn/Biblatex_citation_styles

# Post-processing

Organization:
We want to have post-processing tasks,
which include generating plots from the available plotdata.

Options:
1) could scan the filesystem for plotdata files, and then generate the plots.
But that doesn't tell you how to construct any plots with more than one series on them.
2) Could just define a bunch of scripts, which generate plots.
This is the most flexible approach.
And we probably will do this, at least for things like the brainy-media figure.
But, what about plots that are fairly repetitive, just on different data?
3) Could create plot definition files, kind of like I have before.

Let's get specific. What plots do I have so far to generate?
- The brainy-media figure. This probably belongs in its own script.
- The centerline plots (concentration, and now also electric potential)

So, the question at the moment is, how do I want the centerline plots to work?
Each one is potentially from a different mesh (although not all are).
Do you want series from the same mesh on the same axes?

It sounds like, for the time being, we just need to use scripts.
As we get those developed, some (but not all) may become doit tasks,
and some (but not all) may refactor into data-driven approaches.


_EFFORT_ Function to extract data for 1D plot
- needs to return two 1D arrays (independent and dependent variables)
- could give it coordinate axis, limits, and number of points
- but more generally, could give it two points, and a number of points
- in some cases you may want to put this in a loop, to run multiple lines and plot them on the same set of axes
- could also want different components of a vector on the same set of axes
- need to restrict to points that are actually inside the mesh or on the boundary?

_EFFORT_
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

# Parametric variations
- a given volume fraction can be obtained for different cell and pore sizes, but we can probably just stick with the ones similar to the physical measurements
- the length of the pore, although again we'll probably stick to physical measurements
- the bulk space above and below the pore (far enough away to not affect results)
- mesh refinement study, of course

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
