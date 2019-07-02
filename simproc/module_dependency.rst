
Module Dependency Graph
################################################################################

This is only a listing of the dependencies of these modules on each other,
not on the standard library or external libraries.

Keeping track of this is important in order to help prevent circular dependencies.

- package ``requesthandler``
  - ``filepath``: None
  - ``yaml_manager``: None
  - ``pickle_manager``: None
  - ``nested``: None
  - ``locators``: ``filepath``, ``yaml_manager``
  - ``request``: ``filepath``, ``yaml_manager``, ``locators``, ``nested``
  - ``requestfile``: ``filepath``, ``yaml_manager``, ``locators``, ``request``
  - ``simultaneous``: ``filepath``, ``yaml_manager``, ``request``, ``locators``
  - ``mpi_run``: ``filepath``, ``request``, ``yaml_manager``, ``simultaneous``, ``locators``
  - ``shell``: ``request``, ``yaml_manager``
  - ``debug``: ``request``, ``yaml_manager``, ``shell``
  - ``cleanup``: ``locators``, ``request``, ``yaml_manager``
  - ``customization``: ``filepath``, ``locators``, ``request``, ``yaml_manager``
  - ``commandseq``: ``yaml_manager``, ``pickle_manager``, ``customization``
  - ``templates``: ``yaml_manager``, ``customization``
  - ``generate``: ``yaml_manager``, ``request``, ``customization``
  - ``__init__``: ``filepath``, ``locators``, ``requestfile``, ``customization``, ``shell``, ``templates``, ``cleanup``, ``comparison, ``debug``
  - ``cmdline``: ``*`` (meaning everything listed in ``__init__``)
  - ``__main__``: ``cmdline``

- package ``meshgen``
  - ``gmsh_runner``: ``requesthandler``
  - ``dconv_runner``: ``requesthandler``
  - ``hdf5_conv``: ``requesthandler``, ``dconv_runner``
  - ``onestep``: ``requesthandler``, ``gmsh_runner``, ``dconv_runner``, ``hdf5_conv``

- package ``simulation``
  - ``unitsystem``: None
  - ``meshinfo``: ``requesthandler``
  - ``equationbuilder``: None
  - ``plotseries``: None
  - ``simrequest``: ``requesthandler``, ``meshinfo``, ``equationbuilder``, ``plotseries``
  - ``fickian_homog``: ``requesthandler``, ``meshinfo``, ``equationbuilder``, ``simrequest``

- package ``postproc``
  - ``collection``: ``requesthandler``
  - ``plotting``: ``requesthandler``, ``plotseries``
