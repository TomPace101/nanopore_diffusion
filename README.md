
# Description

Source code and problem description document for FEM solution of diffusion through a nanoporous membrane.

# Setup

__TODO__

# Contacts

- Tom Pace
- Dr. Pete Kekenes-Huskey

# Files and Folders

Otherwise empty directories (eg directories containing only untracked files or subdirectories)
will have a '.keep' file to force git to include the (nearly) empty directory itself.

- README.md: this document
- notes.md: random notes, thoughts, todo lists, and development plans
- description: problem description document
    - fig_svg: inkscape drawings for figures
    - fig_pdf: figures converted to pdf format
- src: code for FEM analysis
    - mesh: code and data for mesh generation using gmsh
        - geo: storage place for gmsh .geo files
        - gmsh_out: storage place for gmsh terminal output
        - msh: storage place for gmsh .msh files
        - xml: storage place for xml files readable by FEniCS
    - params: storage place for parameter sets
    - postproc: post-processing
    - solutions: stored FEM results (one subdirectory per analysis name)
    - solver: code using FEniCS to generate solutions
