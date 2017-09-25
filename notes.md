Current Goal: x=free volume fraction, y=Deff/Dbulk

_ACTION_ look up Fickian diffusion equation in Auriault.
_ACTION_ look at overleaf shared by Pete
_ACTION_ compile TODO lines from all code files

# Figure
- need to compute effective diffusion constant

# Post-processing
How can I generate plots?
https://fenicsproject.org/qa/11876/extract-solution-at-a-set-of-nodes

# Problem Description
continue developing document.

# Mesh
The doit file currently has some bugs. See TODO list in there.

_EVENTUALLY_ post the jinja2 templates and related code to labscripts

## Storage
Is there a logical way to store created msh/xml files for later retrieval if everything is the same?
doit seems to be the way to do this.
It could associate a dictionary of the parameters, including a hash of the relevant template(s), with each stored file.
Or maybe those are separate things: the file depends on the template(s) used to generate it,
but the dictionary of parameters otherwise just selects the particular file.
So the mapping needed is from a set of parameters (probably a named tuple) to a filename.
But if doit is going to check for changes to dependencies, it has to know the mapping to filenames.
So that mapping can be in buildgeom, but only if used as a module, not a script.

## Parametric variations
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
