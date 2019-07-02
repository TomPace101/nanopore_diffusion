"""Defining the data series that go on plots"""

#Standard library
from __future__ import print_function, division #Python 2 compatibility
from copy import deepcopy

#Site packages
import numpy as np

#This package
from ..requesthandler import yaml_manager
from ..requesthandler.commandseq import CommandSequenceRequest

class PlotSeries(object):
  """Data for a single series on a plot

  Attributes:

    - xvals = array of x-values
    - yvals = array of y-values
    - label = legend label
    - metadata = other parameters needed to identify the data series
  """
  _expected_attrs=['xvals','yvals','label','metadata']
  def __init__(self,xvals,yvals,label=None,metadata=None):
    for attrname in self._expected_attrs:
      setattr(self,attrname,locals()[attrname])
  def add_to_axes(self,ax,fmt,newlabel=None,**kwd):
    """Plot this series on the specified axes

    This is a wrapper for ax.plot

    Arguments:

      - ax = matplotlib Axes object
      - fmt = matplotlib format string
      - newlabel = label to use for the legend instead of the series-provided one
      - \**kwd = other keyword arguments for Axes.plot

    Returns:

      - The result of call to ax.plot"""
    if newlabel is None:
      label=getattr(self,'label',None)
      if label is None:
        label=''
    else:
      label=newlabel
    return ax.plot(self.xvals,self.yvals,fmt,label=label,**kwd)
  def __getstate__(self):
    """Used for pickling, and converting to yaml"""
    state={}
    for attrname in self._expected_attrs:
      state[attrname]=getattr(self,attrname)
    return state
  def __setstate__(self,state):
    """Used for unpickling, and loading from yaml"""
    self.__init__(**state)

class DefinePlotSeries(CommandSequenceRequest):
  """A class to load, generate, modify, and save PlotSeries objects using command sequences"""
  def from_dataframe(dfpath,xcol,ycol,outattr,query=None,label=None,metadata=None):
    """Generate a series based on previously-loaded pandas dataframe

    Arguments:
    
      - dfpath = nested path to the DataFrame
      - xcol = name of the DataFrame column to use for the x-values of the series
      - ycol = name of the DataFrame column to use for the y-values of the series
      - outattr = nested path for storing the resulting series
      - query = optional, DataFrame query to use to select only certain rows of the dataframe, as an argument to pd.DataFrame.query
      - label = optional series label
      - metadata = optional series metadata
    
    If you want to save the result to a pickle file, use a subsequent call to ``save_pickle``."""
    #Get the dataframe
    df = self.get_nested(dfpath)
    #Select particular rows if requested
    if query is None or len(query)==0:
      qdf=df
    else:
      qdf=df.query(query)
    #Get the data for the series from the dataframe
    xvals=qdf[xcol]
    yvals=qdf[ycol]
    #Instantiate the series
    series=PlotSeries(xvals=vals,yvals=yvals,label=label,metadata=metadata)
    #Store result
    self.set_nested(outattr,series)
  def from_normalization(inpath,outpath,normpath,label=None,metadata=None):
    """Generate a PlotSeries by normalizing another one

    Arguments:

      - inpath = nested path to the input data series
      - outpath = nested path to the output data series
      - normpath = nested path to the normalization constant
      - label = optional series label
      - metadata = optional series metadata (added to the metadata from the other series)"""
    #Data from other series
    other=self.get_nested(inpath)
    xvals=other.xvals
    orig_y=np.array(other.yvals)
    out_meta=deepcopy(other.metadata)
    #Do the normalization
    factor=self.get_nested(normpath)
    yvals=orig_y/factor
    #Update metadata
    out_meta.update(metadata)
    #Instantiate the series
    series=PlotSeries(xvals=vals,yvals=yvals,label=label,metadata=out_meta)
    #Store result
    self.set_nested(outpath,series)

#Register for loading from yaml
yaml_manager.register_classes([PlotSeries, DefinePlotSeries])
