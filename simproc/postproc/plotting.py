"""Generate plots"""

#Standard library

#Site packages
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd

#This package
from ..requesthandler.commandseq import WithCommandsRequest, make_schema
from ..requesthandler import yaml_manager
from .plotseries import PlotSeries

_FigureRequest_props_schema_yaml="""#FigureRequest
loadfiles:
  type: object
  additionalProperties:
    items: pathlike
rcparams: {type: object}
prepcommands: {type: array}
plotcommands: {type: array}"""


class FigureRequest(WithCommandsRequest):
  """For generating matplotlib figures
  
  User-defined attributes:
  
    - loadfiles = dictionary of pickle files to load (yaml files may be loaded using a prepcommand), usually to load series for ploting
        The dictionary is formatted as {attribute: file path or locator},
        where the attribute to store the loaded data could also be an attribute path.
    - rcparams = dictionary of matplotlib rcParams to be set
    - prepcommands = sequence of commands to execute prior to plotting, such as for defining series (see below)


    prepcommands and plotcommands are both command lists.
    Each command in the list is a pair (cmdname, arguments), where:

      - cmdname = name of the object's method to call, as a string
      - arguments = dictionary of arguments to the method: {argname: value,...}    
  """
  _self_task=True
  _config_attrs=('prepcommands','plotcommands')
  _required_attrs=[]
  _props_schema=make_schema(_FigureRequest_props_schema_yaml)
  def __init__(self,**kwargs):
    #Initialization from base class
    super(FigureRequest, self).__init__(**kwargs)
    #Default command attributes
    if not hasattr(self,'prepcommands'):
      self.prepcommands=[]
    if not hasattr(self,'plotcommands'):
      self.plotcommands=[]
    #Get input files
    self._more_inputfiles=getattr(self,'_more_inputfiles',[]) #Initialize attribute if it doesn't already exist
    self._more_inputfiles+=getattr(self,'loadfiles',[])
    self._more_inputfiles+=self.list_iofiles(self.prepcommands,['filename','infpath'],'_inputfiles')
    self._more_inputfiles+=self.list_iofiles(self.plotcommands,['filename','infpath'],'_inputfiles')
    #Get output files
    self._more_outputfiles=getattr(self,'_more_outputfiles',[]) #Initialize attribute if it doesn't already exist
    self._more_outputfiles+=self.list_iofiles(self.prepcommands,['filename','outfpath'],'_outputfiles')
    self._more_outputfiles+=self.list_iofiles(self.plotcommands,['filename','outfpath'],'_outputfiles')
  def run(self):
    #Final checks and preparatory steps
    self.pre_run()
    #Load data files
    for attrpath,floc in getattr(self,'loadfiles',{}).items():
      self.load_pickle(floc,attrpath)
    #Run prepcommands
    self.process_command_sequence(attrpath='prepcommands',singlefunc=None,positional=False)
    #Set rcParams
    self.set_rcparams(**getattr(self,'rcparams',{}))
    #Set up figures
    ##TODO
    #Set up axes
    ##TODO
    #Add series to axes
    ##TODO
    #Run plotcommands
    self.process_command_sequence(attrpath='plotcommands',singlefunc=None,positional=False)
    #Save figures
    ##TODO
    #Done
    return
  def series_from_dataframe(self,dfpath,xcol,ycol,outattr,query=None,label=None,metadata=None):
    """Generate a series based on previously-loaded pandas dataframe

    Arguments:
    
      - dfpath = nested path to the DataFrame (NOT a file path or locator)
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
    series=PlotSeries(xvals=xvals,yvals=yvals,label=label,metadata=metadata)
    #Store result
    self.set_nested(outattr,series)
  def series_from_normalization(self,inpath,outpath,normpath,label=None,metadata=None):
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
  def set_rcparams(self,**kwargs):
    """To set matplotlib rcParams"""
    for k,v in kwargs.items():
      mpl.rcParams[k]=v
    return
  def figure(self,figattr="fig",**kwargs):
    """To create a new matplotlib.figure instance

    Arguments:

      - figattr = nested path to store generated figure, defaults to "fig"
      - \**kwargs = keyword arguments to pass to plt.figure"""
    fig=plt.figure(**kwargs)
    self.set_nested(figattr,fig)
    return
  def savefig(self,outfpath,figattr="fig",**kwargs):
    """Save the figure to file"""
    fig=self.get_nested(figattr)
    fig.savefig(self.renderstr(outfpath),**kwargs)
  def closefig(self,figattr="fig"):
    """Close the requested figure"""
    fig=self.get_nested(figattr)
    plt.close(fig)
  def axes(self,nrows=1,ncols=1,index=1,figattr="fig",axattr="ax",**kwargs):
    """To create a new matplotlib Axes instance

    Arguments:

      - nrows = number of rows of subplots, defaults to 1
      - ncols = number of columns of subplots, defaults to 1
      - index = index number of subplot, identifying which row and column (see matplotlib docs for explanation), defaults to 1
      - figattr = nested path to figure instance, defaults to "fig"
      - axattr = nested path to resulting Axes instance, defaults to "ax"
      - \**kwargs = keyword arguments to pass to plt.figure.subplots"""
    ax=self.get_nested(figattr).add_subplot(nrows,ncols,index,**kwargs)
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
