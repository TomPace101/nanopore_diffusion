
#Standard library
import os
import os.path as osp
import sys

#Site packages
import numpy as np

#Local
import useful


class PlotSeries(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store the data for a single series on a plot.
  Attributes:
    xvals = array of x-values
    yvals = array of y-values
    label = legend label
    metadata = other parameters needed to identify the data series
  """
  __slots__=['xvals','yvals','label','metadata']
  def add_to_axes(self,ax,fmt,**kwd):
    """Plot this series on the specified axes
    This is a wrapper for ax.plot
    Arguments:
      ax = matplotlib Axes object
      fmt = matplotlib format string
      **kwd = other keyword arguments for Axes.plot
    Returns:
      The result of call to ax.plot"""
    return ax.plot(self.xvals,self.yvals,fmt,label=self.label,**kwd)
