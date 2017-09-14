Current Goal: x=free volume fraction, y=Deff/Dbulk

_ACTION_ set up doit for rebuilding the image pdfs and description pdf (see log 2017-09-07)
_ACTION_ look up Fickian diffusion equation in Auriault.
_ACTION_ look at overleaf shared by Pete
_ACTION_ compile TODO lines from all code files

# Steps
- first, write down the problem
- then, code an outline of the solution
- use doit to break into steps, with output files as the connecting pieces
- borrow from the old code where applicable
- do the standard diffusion equation first, them smoluchowski
- maybe then do a charged solute but a fixed potential, then do the applied charge BCs.
- start with a contiguous domain (a single hole) first,
- then try to do a unit cell

# Problem Description
_ACTION_ create doit file
Then, continue developing document.

# Mesh
_ACTION_ doit automation (make geo, make msh, make xml)

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