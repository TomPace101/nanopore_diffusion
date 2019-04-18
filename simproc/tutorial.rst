
.. command-line usage: python -m doctest  tutorial.rst

Tutorial
################################################################################

Testing
=======

Here are the commands that I use right now to test everything.
From within ``pore``:

- ``python -m simproc --validate``
- ``python -m simproc validation/all_validation.yaml``
- ``doit control=validation/all_validation.yaml``
- ``doit -n 4 control=validation/all_validation.yaml``

Cleanup:

- ``python -m simproc validation/cleanup.yaml``

Requests
========

A Request defines a particular task to be executed.
In particular, a Request instance should be able to do the following three things:

1. Run itself: actually perform the requested action.
   If a request is thought of as a function, this means calling the function.
2. Provide ``doit`` task definitions for its actions.
   The python package ``doit`` provides functionality similar to ``make``,
   but with greater flexibility.
   That is, ``doit`` tracks the state of defined tasks,
   so that only tasks that need to be run will actually be run.
   There are other advantages as well; see http://pydoit.org/ .
3. Validate its own configuration.
   If a request is thought of as a function, this means checking that the input arguments provided are valid.

The process of breaking a simulation up into Requests involves separating the various steps in the process.
Some typical steps are:
- mesh generation
- simulation execution
- post-processing

**TODO**

_TODO_ be sure to describe the ``complete`` method.

.. doctest::
  
  Very basic test of functionality.
  
  >>> import simproc.requesthandler.debug as debug
  >>> req=debug.DummyRequest(name='debug_test',test='abcd1234')
  >>> req.run()
  abcd1234

Locators
========

One of the biggest challenges in a simulation is managing all the data files.
Keeping a consistent file structure makes it easier to find data,
especially when you come back to look at everything after time away.

The ``locators`` module helps to define a consistent set of file and folder names,
and automatically calculate path strings based on this defined structure for a given request.

A locator can be provided for request attribute that is specified as a "Path".

The file structure for the data files starts at a designated "top folder",
which can be indicated with the environment variable ``DATAFOLDER``.
If this environment variable is not defined,
the top folder will be assumed (see module ``locators`` for the assumed location.)
The top folder can also be changed at runtime.

The folder structure within this top folder is determined by the types of locators defined.
For each kind of data file you might want to store, a locator type is created.
Each data file is also assumed to belong to a particular request,
so the locator can use information from the request to calculate the data file path.

**TODO** 

Here is a complete example of the use of locators.

.. doctest::
  
  Setup
  
  >>> #module imports
  >>> from simproc.requesthandler.filepath import Path
  >>> import simproc.requesthandler.locators as locators
  >>> import simproc.requesthandler.debug as debug
  >>> import simproc.requesthandler.yaml_manager as yaml_manager
  >>> #data to be used in multiple examples below
  >>> locators.DATAFOLDER=Path("data") #set the top folder location
  >>> req=debug.DummyRequest(name='debug.alpha.example',test='some_data_here') #an example request
  
  Example 1: string specifiers

  >>> locators.folder_structure.update(TestLocator=['testfolder']) #define a new locator type called 'TestLocator'
  >>> locators.TestLocator.__name__
  'TestLocator'
  >>> tl=locators.TestLocator('myfile.stuff') #Use the TestLocator to locate a file called 'myfile.stuff'
  >>> tl.path(req.name) #Get the path to this file for the example request
  Path('data/testfolder/myfile.stuff')
  >>> #In this case, the request didn't actually contribute to the file path at all
  
  Example 2: request name specifiers

  >>> locators.folder_structure.update(TestLocator=[0,'stuff_files',1]) #Modify TestLocator to use parts of the request name
  >>> tl.path(req.name) #Now where is the test file located?
  Path('data/debug/stuff_files/alpha/myfile.stuff')

  Example 3: going beyond the length of the request name

  >>> locators.folder_structure.update(TestLocator=[0,'stuff_files',1,2,3,4,5]) #This would use up to six parts of a request name
  >>> tl.path(req.name)
  Path('data/debug/stuff_files/alpha/example/myfile.stuff')
  >>> #The non-existent portions of the request name are simply ignored

  Example 4: loading locators from yaml, and writing them to yaml
  
  >>> locators.folder_structure.update(TestLocator=['testing'])
  >>> ys1="!TestLocator test.dat"
  >>> loc=yaml_manager.readstring(ys1)
  >>> loc.path("This string won't appear in the path because of the locator definition")
  Path('data/testing/test.dat')
  >>> ys2=yaml_manager.writestring(loc)
  >>> loc2=yaml_manager.readstring(ys2)
  >>> loc2.path("Again, this string doesn't matter.")
  Path('data/testing/test.dat')
  
  Example 5: defining new locators (or modifying old ones) from yaml
  
  >>> #Define two new locators: InputFile and OutputFile
  >>> #First, let's show that they don't alraedy exist
  >>> locators.InputFile
  Traceback (most recent call last):
    ...
  AttributeError: module 'simproc.requesthandler.locators' has no attribute 'InputFile'
  >>> ys="""!UpdateFolderStructure
  ... InputFile: [input,0]
  ... OutputFile: [output,0]"""
  >>> obj=yaml_manager.readstring(ys)
  >>> #And now they will exist
  >>> loc=locators.InputFile("my_input_file.dat")
  >>> loc.path(req.name)
  Path('data/input/debug/my_input_file.dat')
  >>> loc=locators.OutputFile("my_output_file.dat")
  >>> loc.path(req.name)
  Path('data/output/debug/my_output_file.dat')
  
  Example 6: changing the data folder from within yaml
  
  >>> ys="""!SetDataFolder
  ... datafolder: newdatafolder
  ... resolve: False""" #Normally you would not include this line, but we don't want an absolute path for this example.
  >>> obj=yaml_manager.readstring(ys)
  >>> loc.path(req.name)
  Path('newdatafolder/output/debug/my_output_file.dat')
  
More examples of the use of locators from within yaml files can be found in the validation files.

Miscellany
==========

Just few things not to forget, until I can find a better place for them.

One way to run requests in parallel is to let doit execute the tasks in parallel, with its ``-n`` switch.
Try it for yourself: ``doit -n 4 control=<<requestfile>>``.
There is also now a request subclass that can execute child requests in parallel.

.. doctest::

  Test taking filepath.Path instances round-trip through yaml.

  >>> import simproc.requesthandler.yaml_manager as yaml_manager
  >>> import simproc.requesthandler.filepath as filepath
  >>> p=filepath.Path('/nonexist.txt')
  >>> ys=yaml_manager.writestring(p)
  >>> p2=yaml_manager.readstring(ys)
  >>> p2==p
  True
