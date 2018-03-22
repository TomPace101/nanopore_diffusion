
# Description

Source code and problem description document for FEM solution of diffusion through a nanoporous membrane.

# Required Software

The FEM analysis is conducted using `FEniCS` (https://fenicsproject.org/),
using the python programming language (https://www.python.org/).

Meshes for the FEM analysis are generated using `gmsh` (http://gmsh.info/).

The FEM analysis and the generation of the problem description document
have been automated using a python package known as `doit` (http://pydoit.org/).
The analyses can be run without `doit`, so it is not required.

Many analysis parameters are stored in `yaml` format (https://en.wikipedia.org/wiki/YAML).
A very brief introduction to the syntax can be found at https://learnxinyminutes.com/docs/yaml/,
and its official site is http://yaml.org/.
The python package `pyyaml` (https://github.com/yaml/pyyaml) is used to read these files.

Some output files (such as `gmsh` .geo input files)
are generated using the `jinja2` template engine (http://jinja.pocoo.org/).

Some output data is stored in `pandas` DataFrame objects (http://pandas.pydata.org/).

__TODO__: `matplotlib`

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
    - unitsystem.py: convenience module for converting values to and from model units
    - common.py: functions and classes used by many of the other modules
    - buildgeom.py: code for generating `gmsh` .geo files from input data
    - solver_general.py: functions and classes used by many of the solver modules
    - solver_run.py: top-level module for running solvers
    - solvers: equation-specific solver modules:
        - solver_fickian_unhomog.py: solver for Fick's Law, unhomogenized
        - solver_smol_unhomog.py: solver for Smoluchowski equation, unhomogenized
        - solver_tdpnp_unhomog.py: solver for time-dependent Poisson-Nernst-Planck equation, unhomogenized
        - __TODO__: more to come
    - collect_results.py: generate `pandas` DataFrame from FEM results, and store.
    - plotdata.py: for generating plots
    - postproc.py: code for handling post-processing requests
    - dodo.py: `doit` input file for the FEM analysis
    - paramgen.py: code for generating parameter sets from templates and input data
    - customizations: modules used to customize the equation solvers for specific problems
- data: input and output data from FEM analysis
    - control.yaml: list of analyses to be run with `doit` (see description below)
    - mesh: data for mesh generation using `gmsh`
        - geomdef: geometry defintion files as inputs to mesh generation
        - templates: `jinja2` templates of `gmsh` .geo files
        - geo: storage place for `gmsh` .geo files
        - gmsh_out: storage place for `gmsh` terminal output
        - msh: storage place for `gmsh` .msh files
        - xml: storage place for xml files readable by `FEniCS`
    - paramgen: storage place for templates used in parameter generation
    - params: storage place for parameter sets
        - mesh: storage place for mesh definition parameters
        - model: storage place for model definition parameters
        - paramgen: storage place for definitions of sets of other parameters
        - postproc: storage place for post-processing definition parameters
    - solutions: stored FEM results (one subdirectory per analysis name)
    - postproc: figures generated by post-processing

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

# Usage

__TODO__
- describe control.yaml and how it works with doit

Misc. things to note:
- the distinction between loading objects and actually running them:
  - generally, you could think of it as initialization reads in everything that needs to be read, and running is when output is actually generated
  - also, any quick and simple calculations can be done at initialization, but any time-consuming ones need to wait for running
  - for classes derived from common.ParameterSet, initialization gives all the information needed to create a task in doit, running actually performs the task
  - for classes derived from solver_general.GenericSolver, initialization stops just before asking FEniCS to solve the problem, running solves and generates output
- solver_general defines base classes used by other solver_modules, and solver_run uses those to define tasks and run analyses.
- how files are located (e.g. the `_folders` attribute, etc.)

Discuss how to use unitsystem.py

# Contacts

- Tom Pace
- Dr. Pete Kekenes-Huskey
