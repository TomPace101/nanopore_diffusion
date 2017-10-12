
# Description

Source code and problem description document for FEM solution of diffusion through a nanoporous membrane.

# Required Software

The FEM analysis is conducted using `FEniCS` (https://fenicsproject.org/),
using the python programming language (https://www.python.org/).

Meshes for the FEM analysis are generated using `gmsh` (http://gmsh.info/).

The FEM analysis and the generation of the problem description document
have been automated using a python package known as `doit` (http://pydoit.org/).

Many analysis parameters are stored in `yaml` format (https://en.wikipedia.org/wiki/YAML).
The python package `pyyaml` (https://github.com/yaml/pyyaml) is used to read these files.

Some output files (such as `gmsh` .geo input files)
are generated using the `jinja2` template engine (http://jinja.pocoo.org/).

The problem description document requires the following software:

- `inkscape` (https://inkscape.org/en/)
- `pdflatex` (available with most distributions of `LaTeX`, https://www.latex-project.org/)
- python package `doit`, as described above.

This is a `git` repository (https://git-scm.com/).

# Files and Folders

Otherwise empty directories (eg directories containing only untracked files or subdirectories)
will have a '.keep' file to force `git` to include the directory itself.

- README.md: this document
- notes.md: random notes, thoughts, todo lists, and development plans
- description: problem description document
    - fig_svg: `inkscape` drawings for figures
    - fig_pdf: figures converted to pdf format
- src: code for FEM analysis
    - mesh: code and data for mesh generation using `gmsh`
        - geo: storage place for `gmsh` .geo files
        - gmsh_out: storage place for `gmsh` terminal output
        - msh: storage place for `gmsh` .msh files
        - xml: storage place for xml files readable by `FEniCS`
    - params: storage place for parameter sets
        - control: storage place for run definition parameters
        - scripts: python scripts to generate various parameter definition files
        - mesh: storage place for mesh definition parameters
        - model: storage place for model defintion parameters
    - postproc: post-processing
    - solutions: stored FEM results (one subdirectory per analysis name)
    - solver: code using `FEniCS` to generate solutions

The various sets of analyses that have been conducted using the code are organized
by short, human-readable hashes of longer descriptions of their purposes.
The human readable hashes used here consist of an adjective followed by a noun.
Some readers may be familar with the code names used for releases of the Ubuntu distribution of linux,
which uses this pattern as well.
The human-readable hashes here were generated using the python package `codenamize`
(https://github.com/jjmontesl/codenamize).
The code itself does not use this package, so it is not required
in order to run analyses or generate the problem description document.

# Contacts

- Tom Pace
- Dr. Pete Kekenes-Huskey

