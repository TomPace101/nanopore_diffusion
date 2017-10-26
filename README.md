
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

Some output data is stored in `pandas` DataFrame objects (http://pandas.pydata.org/).

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
    - folderstructure.py: provides information on the folder structure described here to the other modules
    - useful.py: functions and classes used by many of the other modules
    - buildgeom.py: code for generated `gmsh` .geo files from input data
    - solver_general.py: functions and classes used my many of the solver modules
    - equation-specific solver modules:
        - fickian_unhomog.py
    - collect_results.py: generate `pandas` DataFrame from FEM results, and store.
    - plotdata.py: for generating plots
    - dodo.py: `doit` input file for the FEM analysis
    - tasks_mesh.py: `doit` task definitions for mesh generation
    - tasks_solver.py: `doit` task definitions for FEM analysis and data extraction
    - tasks_postproc.py: `doit` task defintions for postprocessing
- data: input and output data from FEM analysis
    - mesh: data for mesh generation using `gmsh`
        - geomdef: geometry defintion files as inputs to mesh generation
        - templates: `jinja2` templates of `gmsh` .geo files
        - geo: storage place for `gmsh` .geo files
        - gmsh_out: storage place for `gmsh` terminal output
        - msh: storage place for `gmsh` .msh files
        - xml: storage place for xml files readable by `FEniCS`
    - params: storage place for parameter sets
        - scripts: python scripts to generate various parameter definition files
        - mesh: storage place for mesh definition parameters
        - model: storage place for model definition parameters
    - solutions: stored FEM results (one subdirectory per analysis name)

The various sets of analyses that have been conducted using the code are organized
by short, human-readable hashes of longer analysis descriptions.
The human readable hashes used here consist of an adjective followed by a noun.
Some readers may be familiar with the code names used for releases of the Ubuntu distribution of linux,
which uses this pattern as well.
The human-readable hashes here were generated using the python package `codenamize`
(https://github.com/jjmontesl/codenamize).
The code itself does not use this package, so it is not required
in order to run analyses or generate the problem description document.
And in some cases, the descriptions have been revised after the hash was generated.
The mapping between the current analysis descriptions and the
hashes is stored in `src/params/hashes.yaml`.

# Contacts

- Tom Pace
- Dr. Pete Kekenes-Huskey

