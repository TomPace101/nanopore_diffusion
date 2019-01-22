
Module Dependency Graph
################################################################################

This is only a listing of the dependencies of these modules on each other,
not on the standard library or external libraries.

Keeping track of this is important in order to help prevent circular dependencies.

- package ``requesthandler``
  - ``filepath``: None
  - ``yaml_manager``: None
  - ``locators``: ``filepath``, ``yaml_manager``
  - ``request``: ``filepath``, ``yaml_manager``, ``locators``
  - ``requestfile``: ``filepath``, ``yaml_manager``, ``locators``, ``request``
  - ``shell``: ``request``, ``yaml_manager``
  - ``debug``: ``request``, ``yaml_manager``, ``shell``
  - ``cleanup``: ``locators``, ``request``, ``yaml_manager``
  - ``customization``: ``filepath``, ``locators``, ``request``
  - ``templates``: ``yaml_manager``, ``customization``
  - ``generate``: ``yaml_manager``, ``request``, ``customization``
  - ``__init__``: ``filepath``, ``locators``, ``requestfile``, ``customization``, ``shell``, ``templates``, ``cleanup``, ``comparison, ``debug``
  - ``cmdline``: ``*`` (meaning everything listed in ``__init__``)
  - ``__main__``: ``cmdline``

- package ``meshgen``
  - ``gmsh_runner``: ``requesthandler``
  - ``dconv_runner``: ``requesthandler``

- package ``simulation``
  - `` ``: