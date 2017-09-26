Current Goal: x=free volume fraction, y=Deff/Dbulk

_ACTION_ look up Fickian diffusion equation in Auriault.
_ACTION_ look at overleaf shared by Pete
_ACTION_ compile TODO lines from all code files

# Figure
- need to compute effective diffusion constant

# Post-processing
_ACTION_ Talk with Ryan about immersed surfaces for post-processing

How can I generate plots?
https://fenicsproject.org/qa/11876/extract-solution-at-a-set-of-nodes

# Problem Description
continue developing document.

# Mesh
_ACTION_ Look up "field" in gmsh tutorials to try to resolve issue with centerline.

The doit file currently has some bugs. See TODO list in there.

_EVENTUALLY_ post the jinja2 templates and related code to labscripts

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


# Parametric variations
- a given volume fraction can be obtained for different cell and pore sizes, but we can probably just stick with the ones similar to the physical measurements
- the length of the pore, although again we'll probably stick to physical measurements
- the bulk space above and below the pore (far enough away to not affect results)
- mesh refinement study, of course

# Steps
- first, write down the problem
- then, code an outline of the solution
- use doit to break into steps, with output files as the connecting pieces
- borrow from the old code where applicable
- do the standard diffusion equation first, then smoluchowski
- maybe then do a charged solute but a fixed potential, then do the applied charge BCs.
- start with a contiguous domain (a single hole) first,
- then try to do a unit cell
