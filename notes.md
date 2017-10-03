
# Figure
Current Goal: x=free volume fraction, y=Deff/Dbulk
- need to compute effective diffusion constant

# Simulation
**PRIORITY** calculate flux integrals, and then effective diffusion constants
_THEN_ we need to run for face-centered geometry as well.

Where should the boundary condition parameters be stored?
They should be separate from the mesh parameters,
because we could run the same mesh with many different boundary conditions.
Or maybe they should be the same, and I just programmatically generate multiple yaml documents.

I think some sort of master table is going to be the only way, probably.

_THEN_ work on homogenized Fickian diffusion equation (see discussion in Auriault)

# Separate solution and post-processing

We don't want to re-run the whole model just to change a plot.
So we need to separate out the plot data.
It can be extracted when you generate the solution.
But we need a way to store it an organize it,
so it can get to the plot that needs it.

Where should we store the effective diffusion constants?
This brings back the idea of a need for a master table of results.


# Post-processing
_EFFORT_ Function to extract data for 1D plot
- needs to return two 1D arrays (independent and dependent variables)
- could give it coordinate axis, limits, and number of points
- but more generally, could give it two points, and a number of points
- in some cases you may want to put this in a loop, to run multiple lines and plot them on the same set of axes
- could also want different components of a vector on the same set of axes
- need to restrict to points that are actually inside the mesh or on the boundary?

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

# Problem Description

continue developing document.
It has its own TODO list, but there is more than that.

For the part where we have figures of the geometry,
it would be nice if these could be auto-generated from the lattice yaml file.
That way, if I add internal surfaces, the drawings could auto-update.

Of course, right now the yaml file doesn't contain the physical dimensions.


# Mesh
_ACTION_ Look up "field" in gmsh tutorials to try to resolve issue with centerline. (tutorial 10, and the manual discussion on controlling mesh size)
_ACTION_ check for compatibility of gmsh version (eg Ruled Surface vs Surface)
_ACTION_ add validation of geometric inputs

_EVENTUALLY_ post the jinja2 templates and related code to labscripts

# Parametric variations
- a given volume fraction can be obtained for different cell and pore sizes, but we can probably just stick with the ones similar to the physical measurements
- the length of the pore, although again we'll probably stick to physical measurements
- the bulk space above and below the pore (far enough away to not affect results)
- mesh refinement study, of course

--------------------------------------------------------------------------------
# Non-active items only below this point

# Steps
- first, write down the problem
- then, code an outline of the solution
- use doit to break into steps, with output files as the connecting pieces
- borrow from the old code where applicable
- do the standard diffusion equation first, then smoluchowski
- maybe then do a charged solute but a fixed potential, then do the applied charge BCs.
- start with a contiguous domain (a single hole) first,
- then try to do a unit cell

# Storage
Is there a logical way to store created msh/xml files for later retrieval if everything is the same?
doit seems to be the way to do this.
It could associate a dictionary of the parameters, including a hash of the relevant template(s), with each stored file.
Or maybe those are separate things: the file depends on the template(s) used to generate it,
but the dictionary of parameters otherwise just selects the particular file.
So the mapping needed is from a set of parameters (probably a named tuple) to a filename.
But if doit is going to check for changes to dependencies, it has to know the mapping to filenames.
So that mapping can be in buildgeom, but only if used as a module, not a script.

This same storage structure can be used for solutions, although clearly the parameters for the mesh and the parameters for the solver are different.
What we need is a consistent stem name, reflecting the parameters used to generate a file.
For example, the .geo files already contain the lattice, but they need other parameters as well.
Obviously, we don't want a huge filename that contains every parameter. (Been there, done that.)
So, then, we need a mapping.
Which way? Or both ways?
Both ways would be preferable, although that can be simulated by means of a ticket file, or a separate mapping.
That is, the mapping does params->basename and either through file lookup or another mapping you can go basename->params.
Options for params->basename:
- sequential numbering (length likely ok, as long as you allow enough digits, but must manually look for matches to avoid duplication)
- hash (probably too long, risk of collision if shorter, very un-human readable)
- abbreviations (probably too long, not quite as future-proof as a hash)
Actually, there are human-readable hashes:
https://gist.github.com/raineorshine/8d67049c0aaaa082614e417660462fda
Options for basename->params:
- keep a master dictionary (perhaps a dataframe)
- ticket files
I think I prefer the former.
In fact, the way this will really work is that when you're going either way you'll want to use a master lookup.
If you're going params->basename, and you don't find it, you create a new entry.
A dataframe might make a lot of sense in this case, although that would potentially lose the ability to nest the parameters.

Of course, the alternative to this is to assume you're only ever doing one run.
And that should certainly be an option: I just want one analysis, so don't do any name storage.
If you want to do it that way, you just provide it with a basename.

Or, if you don't provide a basename, one will be auto-generated using a hash of your other input parameters.
https://github.com/jjmontesl/codenamize

Note that the same master table can be used as a store of single-parameter results, such as diffusion constants.

The final resolution of this was that the basename is provided in the control.yaml file.
It is not auto-generated.

# Doit
Given the discussion above in Storage, how do we implement this in doit?
We want to be able to give it just a single parameter dictionary, including basename, and have it run that.
Sounds like a job for a yaml input file.
Potentially the top-level parameter file just points to mesh and solver parameter input files.

If I want the basename to be auto-generated, I just leave it blank.

What if I want a slew of runs, such as parametric variations?
Yaml files allow multiple "documents", separated by three dashes.
(In python, use yaml.load_all, which returns a generator, which you probably want inside a list comprehension)

What if I want to auto-generate the parameters in a for-loop?
Then write code to generate all the necessary dictionaries,
then save that out to yaml.

Single parameter file, or separate mesh file from solver file?
Single means more duplication, but less complexity.
Single file it is, then.
It is called "control.yaml"
I made it a symbolic link to another yaml file in the params folder.
Each document in control must have its own basename.

What about this master lookup table, then?
Do we really need that, now?
The control yaml specifies the parameters and the basename,
or the basename is consistently (ie repeatably) calculated if not provided.
All the master table helps do is avoid parametric duplicates with different basenames,
and hash collisions. That could be important, though.

Note that for now, the basenames are not calculated.
And I'm not seeing a good reason for them to be at the moment.

