Current Goal: x=free volume fraction, y=Deff/Dbulk

_ACTION_ set up doit for rebuilding the image pdfs and description pdf (see log 2017-09-07)
_ACTION_ look up Fickian diffusion equation in Auriault.
_ACTION_ look at overleaf shared by Pete
_ACTION_ compile TODO lines from all code files

# Problem Description
_ACTION_ put figure svg->pdf step in doit file
Then, continue developing document.

# Mesh
_ACTION_ doit automation (make geo, make msh, make xml)

## Storage
Is there a logical way to store created msh/xml files for later retrieval if everything is the same?
doit seems to be the way to do this.
It could associate a dictionary of the parameters, including a hash of the relevant template(s), with each stored file.
Or maybe those are separate things: the file depends on the template(s) used to generate it,
but the dictionary of parameters otherwise just selects the particular file.
So the mapping needed is from a set of parameters (probably a named tuple) to a filename.
But if doit is going to check for changes to dependencies, it has to know the mapping to filenames.
So that mapping can be in buildgeom, but only if used as a module, not a script.

There's a larger issue here of the geo and msh filenames being part of the parameter dictionary.
Should that be a separate thing, somehow?

## XML
http://mypages.iit.edu/~asriva13/?page_id=586
see also smolhomog code, which I think does it too.

## Parametric variations
- a given volume fraction can be obtained for different cell and pore sizes, but we can probably just stick with the ones similar to the physical measurements
- the length of the pore, although again we'll probably stick to physical measurements
- the bulk space above and below the pore (far enough away to not affect results)
- mesh refinement study, of course

# Boundary Conditions

# Diffusion Equation

# Steps
- first, write down the problem
- then, code an outline of the solution
- use doit to break into steps, with output files as the connecting pieces
- borrow from the old code where applicable
- do the standard diffusion equation first, then smoluchowski
- maybe then do a charged solute but a fixed potential, then do the applied charge BCs.
- start with a contiguous domain (a single hole) first,
- then try to do a unit cell
