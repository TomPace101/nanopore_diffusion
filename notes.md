

# Code/Misc

_TODO_ complete the "test" analysis, even without middle surface

_TODO_ actually generate the plots we extracted data for (concentration vs centerline)

_TODO_ shell command files for generating .msh and .xml files, which doit then calls?

_TODO_ I no longer like the hashes I'm getting from codenamize
We still need something like that, though.
I experimented in the notebook with some other stuff.

_TODO_ wiki page on flux integration over internal boundary
Should probably document how to do external boundary as well, for comparison.

_TODO_ syncing/mirroring
- /u1 (could use unison)
- shared (could use unison)
- holly and/or dlx (has to be rsync, unless we can build unison on the server)

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

_ACTION_ look at Pete's code provided by Bin
_ACTION_ get boltzmann constant into module?

## Homogenized Fickian

_ACTION_ see discussion in Auriault

# Problem Description

_TODO_ bibliography, with biblatex
https://www.sharelatex.com/learn/Bibliography_management_in_LaTeX
bibtex is older:
https://www.sharelatex.com/learn/Bibliography_management_with_bibtex
http://www.bibtex.org/Using/

_TODO_ Review writeup on expected result.
I didn't include the part about how the area needs to be piecewise constant.
Should this be in there?

_DECISION_ maybe the result figures and discussion should be a separate document.
That way, the problem description could be generated before the analyses are run.

You could even set it up so that the document that requires completed runs
also pulls in the problem description itself,
so in the end you do get a document that has everything.
This would require splitting the description up as well:
- the contents that are part of both documents
- the part that pulls that in as a standalone document

It has its own _TODO_ list.

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

# Post-processing
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

_EFFORT_ solve with flux as another primary variable

_EFFORT_ doctests? some other kind of test?
Maybe just create scripts to do the basic validation on buildgeom and the solvers.

_EFFORT_ mesh refinement study
_EFFORT_ study of required H value


--------------------------------------------------------------------------------
# Non-active items only below this point
