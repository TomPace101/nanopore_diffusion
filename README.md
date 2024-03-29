
# Description

Source code and problem description document for FEM solution of diffusion problems at the nanoscale.

The principal code in this repository is now contained in the `simproc` python module.

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
The python package `ruamel.yaml` (https://pypi.org/project/ruamel.yaml/) is used to read these files.

Some output files (such as `gmsh` .geo input files)
are generated using the `jinja2` template engine (http://jinja.pocoo.org/).

Some output data is stored in `pandas` DataFrame objects (http://pandas.pydata.org/).

`matplotlib` (https://matplotlib.org/) is used for generation of plots.

`scipy` (https://www.scipy.org/) is used only for unit conversions and physical constants.

Documentation of the code is set up to be compiled with Sphinx: http://www.sphinx-doc.org/en/master/.
__TODO__ what other software does sphinx require?
  - latex (listed below)
  - pandoc? (I think this was only if you used markdown)
  
The sphinx documentation is built by running the command `make latexpdf` from within the `sphinx` directory.

__TODO: OLD__
The problem description document requires the following software:

- `inkscape` (https://inkscape.org/en/)
- `pdflatex` (available with most distributions of `LaTeX`, https://www.latex-project.org/)
- python package `doit`, as described above.

The code is hosted in a `git` repository (https://git-scm.com/).

__TODO__ explain the singularity recipe, and provide shub link to the image, if possible.

# Files and Folders

Otherwise empty directories (eg directories containing only untracked files or subdirectories)
will have a '.keep' file to force `git` to include the directory itself.

- README.md: this document
- setup.py: file to facilitate installation using setuptools
- pseudocode: pseudocode for exploring new features
- simproc: code for FEM analysis, see its own documentation
- singularity: recipes for singularity images containing the software dependencies
- sphinx: sphinx documentation folder
- validation: validation data for the simproc module
- old_notes.md and new_notes.md: random notes, thoughts, todo lists, and development plans

The following directories are mostly no longer used
- old-sphinx: sphinx documentation for older code
- description: problem description document
    - fig_svg: `inkscape` drawings for figures
    - fig_pdf: figures converted to pdf format
- src: OLD code for FEM analysis
    - dependencies_test.py: file to test that all required software dependencies are satisfied
    - folderstructure.py: provides information on the folder structure described here to the other modules
    - unitsystem.py: convenience module for converting values to and from model units
    - common.py: functions and classes used by many of the other modules
    - buildgeom.py: code for generating `gmsh` .geo files from input data
    - simulator_general.py: functions and classes used by many of the simulator modules
    - simulator_run.py: top-level module for running simulators
    - simulators: equation-specific simulator modules:
        - simulator_fickian_unhomog.py: simulator for Fick's Law, unhomogenized
        - simulator_smol_unhomog.py: simulator for Smoluchowski equation, unhomogenized
        - simulator_tdpnp_unhomog.py: simulator for time-dependent Poisson-Nernst-Planck equation, unhomogenized
        - __TODO__: more to come
    - collect_results.py: generate `pandas` DataFrame from FEM results, and store.
    - plotdata.py: for generating plots
    - postproc.py: code for handling post-processing requests
    - dodo.py: `doit` input file for the FEM analysis
    - paramgen.py: code for generating parameter sets from templates and input data
    - customizations: modules used to customize the equation simulators for specific problems
- data: OLD input and output data from FEM analysis
    - control.yaml: list of analyses to be run with `doit` (see description below)
    - mesh: data for mesh generation using `gmsh`
        - geomdef: geometry defintion files as inputs to mesh generation
        - templates: `jinja2` templates of `gmsh` .geo files
        - geo: storage place for `gmsh` .geo files
        - gmsh_out: storage place for `gmsh` terminal output
        - msh: storage place for `gmsh` .msh files
        - xml: storage place for xml files readable by `FEniCS`
    - paramgen_tmpl: storage place for templates used in parameter generation
    - params: storage place for parameter sets
        - mesh: storage place for mesh definition parameters
        - model: storage place for model definition parameters
        - paramgen: storage place for definitions of sets of other parameters
        - postproc: storage place for post-processing definition parameters
    - solutions: stored FEM results (one subdirectory per analysis name)
    - postproc: figures generated by post-processing

If defined, the environment variable `SRCLOC` will be used to locate the `src` folder.
Similarly, if the environment variable `DATALOC` is defined, it will be used to locate the `data` folder.

The various sets of analyses in the old `data` folder are organized
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
  - for classes derived from simulator_general.GenericSimulator, initialization stops just before asking FEniCS to solve the problem, running solves and generates output
- simulator_general defines base classes used by other simulator modules, and simulator_run uses those to define tasks and run analyses.
- how files are located (e.g. the `_folders` attribute, etc.)

Discuss how to use unitsystem.py

# Contacts

- Tom Pace
- Dr. Pete Kekenes-Huskey
