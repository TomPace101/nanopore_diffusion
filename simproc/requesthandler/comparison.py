"""For comparing output to the expected output, to validate the code

Setup data for test

>>> files=['/tmp/a.txt', '/tmp/b.txt', '/tmp/c.txt']
>>> files_data=['hello world\\n']*2 + ['goodbye!\\n']
>>> for fpath,data in zip(files,files_data):
...   with open(fpath,'w') as fh:
...     fh.write(data)
...
12
12
9

Example with files that do match

>>> req1=FileComparisonRequest(expected=files[0],received=files[1])
>>> req1.run()
'Files match: /tmp/a.txt and /tmp/b.txt'
>>> req1b=FileSizeComparisonRequest(expected=files[0],received=files[1])
>>> req1b.run()
'File sizes match: /tmp/a.txt and /tmp/b.txt'

Example with files that don't match

>>> req2=FileComparisonRequest(expected=files[1],received=files[2])
>>> req2.run()
Traceback (most recent call last):
  ...
AssertionError: Found unexpected difference in files /tmp/b.txt and /tmp/c.txt
>>> req2b=FileSizeComparisonRequest(expected=files[1],received=files[2])
>>> req2b.run()
Traceback (most recent call last):
  ...
AssertionError: Found unexpected difference in file sizes: /tmp/b.txt has size 12, /tmp/c.txt has size 9, range (0,0) gives limits (12,12)
>>> req2c=FileSizeComparisonRequest(expected=files[1],received=files[2],range=(-3,5))
>>> req2c.run()
'File sizes match: /tmp/b.txt and /tmp/c.txt'

Clean up

>>> for fpath in files:
...   os.remove(fpath)
...

"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import filecmp

#Site packages

#This package
from . import request
from . import yaml_manager

_FileComparisonRequest_props_schema_yaml="""#FileComparisonRequest
expected: {type: pathlike}
received: {type: pathlike}"""

class FileComparisonRequest(request.Request):
  """Request to compare the contents of two files, and issue error if they don't match
  
  User-Provided Attributes:
  
    - expected: path to the file containing the expected output
    - received: path to the file containing the output produced by the code
    
  When run, returns a comparison report as a string if files match.
  Raises AssertionError if they do not."""
  _props_schema=request.make_schema(_FileComparisonRequest_props_schema_yaml)
  _required_attrs=['expected','received']
  _inputfile_attrs=['expected','received']
  _self_task=True
  def run(self):
    args=(str(self.expected),str(self.received))
    ans=filecmp.cmp(*args,shallow=False)
    assert ans, "Found unexpected difference in files %s and %s"%args
    return "Files match: %s and %s"%args

_FileComparisonListRequest_props_schema_yaml="""#FileComparisonListRequest
pairs:
  type: array
  items:
    type: array
    minItems: 2
    maxItems: 2
    items: {type: pathlike}
_children: {type: array}"""

class FileComparisonListRequest(request.Request):
  """Perform file comparison on multiple pairs of files, and issue error if any pair doesn't match
  
  User-Provided Attributes:
  
    - pairs: list of pairs (expected, received), each a path

  Calculated Attributes:
  
    - _children: A list storing all child requests
    
  Collects the reports or errors of the child requests,
  and raises the errors of all children, if any.
  If no errors, returns the reports of all the child requests."""
  _props_schema=request.make_schema(_FileComparisonListRequest_props_schema_yaml)
  _required_attrs=['pairs']
  _child_seq_attrs=['_children']
  _self_task=False #This request generates doit tasks from its children, not itself
  def __init__(self,**kwargs):
    #Initialization from base class
    super(FileComparisonListRequest, self).__init__(**kwargs)
    #Each listed pair defines a FileComparisonRequest
    self._children=[FileComparisonRequest(expected=p[0],received=p[1]) for p in self.pairs]
  def run(self):
    """Run all child requests and capture their results"""
    #Final checks and preparatory steps
    self.pre_run()
    #Storage for reports and errors
    reportlist=[]
    errlist=[]
    #We only need to call run() on the immediate children.
    #Children with their own children will do the same.
    for req in self.all_children():
      try:
        report=req.run()
        reportlist.append(report)
      except AssertionError as einst:
        errlist.append(str(einst))
    #Raise any errors
    assert len(errlist)==0, "\n".join(errlist)
    return "\n".join(reportlist)

_FileSizeComparisonRequest_props_schema_yaml="""#FileSizeComparisonRequest
expected: {type: pathlike}
received: {type: pathlike}
range:
  anyOf:
    - {type: number}
    - {type: array}
"""

class FileSizeComparisonRequest(request.Request):
  """Request to compare the sizes of two files, within a given range
  
  User-Provided Attributes:
  
    - expected: path to the file containing the expected output
    - received: path to the file containing the output produced by the code
    - range: number or pair specifying the range of variation allowed, defaults to zero
        As a pair, the range is specified as (``lower``,``upper``) or (``upper``, ``lower``):
        The smaller of the two values is taken as ``lower``, and the larger one as ``upper``.
        If a number, then both ``lower`` and ``upper`` are set to this number.
    
  To pass, the received file size must be within the range:
    (expected size + ``lower``) to (expected size + ``upper``), inclusive of the endpoints
    
  Note that this means you probably want ``lower`` to be zero or negative.

  When run, returns a comparison report as a string if the check passes.
  Raises AssertionError if not."""
  _props_schema=request.make_schema(_FileSizeComparisonRequest_props_schema_yaml)
  _required_attrs=['expected','received']
  _inputfile_attrs=['expected','received']
  _self_task=True
  err_tmpl="Found unexpected difference in file sizes: %s has size %d, %s has size %d, range (%d,%d) gives limits (%d,%d)"
  def run(self):
    #Get the range
    rg=getattr(self,'range',0)
    if hasattr(rg,'__len__'):
      upper,lower=rg
      if upper<lower:
        upper,lower=(lower,upper)
    else:
      upper=rg
      lower=rg
    #Read the size of both files
    exp_size=os.stat(str(self.expected)).st_size
    rcv_size=os.stat(str(self.received)).st_size
    maxsize=exp_size+upper
    minsize=exp_size+lower
    valtup=(str(self.expected),exp_size,str(self.received),rcv_size,upper,lower,minsize,maxsize)
    assert rcv_size >= minsize and rcv_size <= maxsize, self.err_tmpl%valtup
    return "File sizes match: %s and %s"%(str(self.expected),str(self.received))

#Register for loading from yaml
yaml_manager.register_classes([FileComparisonRequest, FileComparisonListRequest, FileSizeComparisonRequest])
