
.. command-line usage: python -m doctest  requesthandler/tutorial.rst

Tutorial
################################################################################

**TODO** 

More on Locators
================

Here is a complete example of the use of locators.

.. doctest::
  
  Setup
  
  >>> from requesthandler.filepath import Path
  >>> import requesthandler.locators as locators
  >>> import requesthandler.debug as debug
  >>> locators.DataFile.parentpath=Path("data")
  >>> req=debug.DummyRequest(name='debug.alpha.example',test='some_data_here')
  
  Example 1: string specifiers

  >>> locators.folder_structure.update(TestLocator=['testfolder'])
  >>> locators.TestLocator.__name__
  'TestLocator'
  >>> tl=locators.TestLocator('myfile.stuff')
  >>> tl.path(req)
  Path('data/testfolder/myfile.stuff')
  
  Example 2: request name specifiers

  >>> locators.folder_structure.update(TestLocator=[0,'stuff_files',1])
  >>> tl.path(req)
  Path('data/debug/stuff_files/alpha/myfile.stuff')

  Example 3: going beyond the length of the request name

  >>> locators.folder_structure.update(TestLocator=[0,'stuff_files',1,2,3,4,5])
  >>> tl.path(req)
  Path('data/debug/stuff_files/alpha/example/myfile.stuff')
