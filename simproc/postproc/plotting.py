"""Generate plots"""

#Standard library

#Site packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#This package
from ..requesthandler.commandseq import CommandSequenceRequest, make_schema
from ..requesthandler import yaml_manager
from .plotseries import PlotSeries

class FigureRequest():
  """For generating matplotlib figures"""
  def figure(self,figattr="fig",**kwargs):
    """To create a new matplotlib.figure instance

    Arguments:

      - figattr = nested path to store generated figure, defaults to "fig"
      - \**kwargs = keyword arguments to pass to plt.figure"""
    fig=plt.figure(**kwargs)
    self.set_nested(figattr,fig)
    return
  def axes(self,nrows=1,ncols=1,index=1,figattr="fig",axattr="ax",**kwargs):
    """To create a new matplotlib Axes instance

    Arguments:

      - nrows = number of rows of subplots, defaults to 1
      - ncols = number of columns of subplots, defaults to 1
      - index = index number of subplot, identifying which row and column (see matplotlib docs for explanation), defaults to 1
      - figattr = nested path to figure instance, defaults to "fig"
      - axattr = nested path to resulting Axes instance, defaults to "ax"
      - \**kwargs = keyword arguments to pass to plt.figure.subplots"""
    ax=self.get_nested(figattr).subplots(nrows,ncols,index,**kwargs)
    self.set_nested(axattr,ax)
    return
  def iteraxes(self,axlist):
    """Return an iterator over the axes"""
    if axlist is None:
      axlist = ["ax"]
    for axpath in axlist:
      yield self.get_nested(axpath)
  def axmethod(self,method,axlist=None,**kwargs):
    """Call the given matplotlib Axes method on the list of axes, with the provided arguments

    Arguments:

      - method = name of Axes method to call, as string
      - axlist = list of nested attribute paths to the axes on which to apply this command, defaults to ["ax"]
      - \**kwargs = keyword arguments to the method"""
    for ax in self.iteraxes(axlist):
      f=getattr(ax,method)
      if kwargs is None:
        kwd={}
      else:
        kwd=kwargs
      f(**kwargs)
    return
  def add_series(self,seriespath,axlist=None,fmt=None,newlabel=None,**kwargs):
    """Add a single series to the listed axes

    Arguments:

      - seriespath = nested path to the series instance to add to the axes
      - axlist = list of nested paths to the Axes instances, defaults to ["ax"]
      - fmt = matplotlib format string
      - newlabel = label to use for the legend instead of the series-provided one
      - \**kwargs = other keyword arguments for Axes.plot"""
    series=self.get_nested(seriespath)
    for ax in self.iteraxes(axlist):
      series.add_to_axes(ax,fmt,newlabel,**kwargs)
    return
  def add_multi_series(self,serieslist,axlist=None,fmtlist=None,labellist=None):
    """Add multiple series to the listed axes

    Arguments:

      - serieslist = list of nested paths to series instances to add to the axes
      - axlist = list of nested paths to Axes instances, defaults to ["ax"]
      - fmtlist = list of format specifier strings, defaults to ``None`` for each series
      - labellist = list of alternative labels for each """
    if fmtlist is None:
      fmtlist=[None]*len(serieslist)
    assert len(serieslist)==len(fmtlist), "Got %d formats for %d series"%(len(fmtlist),len(serieslist))
    if labellist is None:
      labellist=[None]*len(serieslist)
    assert len(serieslist)==len(labellist), "Got %d labels for %d series"%(len(labellist),len(serieslist))
    for i,seriespath in enumerate(serieslist):
      self.add_series(seriespath,axlist,fmtlist[i],labellist[i],**kwargs)
    return
  def hline(self,valpath,axlist=None,**kwargs):
    """Add a horizontal line to the axes

    Arguments:

      - valpath = nested path to the y-value for the horizontal line
      - \**kwargs = keyword arguments for ax.axhline"""
    yval=self.get_nested(valpath)
    for ax in self.iteraxes(axlist):
      ax.axhline(yval,**kwargs)
    return
  def vline(self,valpath,axlist=None,**kwargs):
    """Add a horizontal line to the axes

    Arguments:

      - valpath = nested path to the x-value for the vertical line
      - \**kwargs = keyword arguments for ax.axvline"""
    xval=self.get_nested(valpath)
    for ax in self.iteraxes(axlist):
      ax.axvline(xval,**kwargs)
    return
  def equalityline(self,serpath,axlist=None,label=None,metadata=None):
    """Create a series for a 1:1 line on the requested axes, based on the axis limits

    Because it's based on the axis limits, you should set those before calling this.

    Arguments:

      - serpath = nested attribute path to store the series at
      - axlist = list of nested paths to Axes instances, defaults to ["ax"]
      - label = label to use for the series
      - metadata = metadata to use for the series"""
    for ax in self.iteraxes(axlist):
      xmin,xmax=ax.get_xlim()
      ymin,ymax=ax.get_ylim()
      vals=(min(xmin,ymin), max(xmax,ymax))
      ser=PlotSeries(xvals=vals,yvals=vals,label=label,metadata=metadata)
      self.set_nested(serpath,ser)
    return

#Register for loading from yaml
yaml_manager.register_classes([FigureRequest])
