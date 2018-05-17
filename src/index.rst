.. Nanoscale Diffusion documentation master file, created by
   sphinx-quickstart on Mon May 14 16:39:45 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Nanoscale Diffusion's documentation!
===============================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

General Overview
================

This collection of scripts makes extensive use of ``yaml`` format (https://en.wikipedia.org/wiki/YAML) for input files,
and even some output files as well.
A very brief introduction to the syntax can be found at https://learnxinyminutes.com/docs/yaml/,
and its official site is http://yaml.org/.

This collection of scripts can perform the following steps:

#. Generate meshes for FEM simulation: MeshGeneration_
#. Run simulations in ``FEniCS`` and extract selected data: Simulation_
#. Post-process the simulation results to generate graphics: PostProc_
#. Generate input parameter files for the steps above from simpler input files: ParamGen_

There are other additional MiscModules_ as well.

**TODO**: list required software, as in README

**TODO**: describe expected folder structure, as in README

.. _MeshGeneration:
Mesh Generation
===============

This process creates a ``.geo`` file, which ``gmsh`` then processes into ``.msh`` format,
which is then converted to ``.xml`` usable by ``FEniCS`` with ``dolfin-convert``.

buildgeom.py
------------

.. automodule:: buildgeom
   :members:

geom_mk_msh.py
--------------

.. automodule:: geom_mk_msh
   :members:

geom_mk_xml.py
--------------

.. automodule:: geom_mk_xml
   :members:

.. _Simulation:
Simulations
===========

Each type of simulation that can be performed has its own simulator module.
These modules all use clases defined in simulator_general_.
The correct simulator module can be run by simulator_run_.

Individual simulations can also customize the behavior of the simulator by requesting modules from ``customizations``.

**TODO** explain how customization works

.. _simulator_general:
simulator_general.py
--------------------

.. automodule:: simulator_general
   :members:

.. _simulator_run:
simulator_run.py
----------------

.. automodule:: simulator_run
   :members:

fickian_unhomog.py
------------------

.. automodule:: fickian_unhomog
   :members:

smol_unhomog.py
---------------

.. automodule:: smol_unhomog
   :members:

tdpnp_unhomog.py
----------------

.. automodule:: tdpnp_unhomog
   :members:

.. _PostProc:
Post-Processing
===============

This step is run by postproc.py_, which can call collect_results_ and plotdata_.

.. _postproc.py:
postproc.py
-----------

.. automodule:: postproc
   :members:

.. _collect_results:
collect_results.py
------------------

.. automodule:: collect_results
   :members:

.. _plotdata:
plotdata.py
-----------

.. automodule:: plotdata
   :members:

.. _ParamGen:
Parameter Generation
====================

paramgen.py
-----------

.. automodule:: paramgen
   :members:

.. _MiscModules:
Miscellaneous Modules
=====================

common.py
---------

.. automodule:: common
   :members:

folderstructure.py
------------------

.. automodule:: folderstructure
   :members:

dependencies_test.py
--------------------

.. automodule:: dependencies_test
   :members:

unitsystem.py
-------------

.. automodule:: unitsystem
   :members:

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
