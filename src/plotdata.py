
#Standard library
import os
import os.path as osp

#Site packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Local
from folderstructure import *
import useful
import solver_general

class PlotSeries(useful.ParameterSet):
  """Data for a single series on a plot
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

class PlotFigure(useful.ParameterSet):
  """Data for a single matplotlib figure
  This is for a plot with a single set of axes.
  Attributes:
    To be read in from yaml file:
      figsize = pair of numbers representing figure size, in inches: (Width, Height)
      filename = name of the output file to be created, as string
      prepfunctions = sequence of method calls used to generate additional data, etc.
        The available method names usually start with 'prep_'
      plotfunctions = sequence of method calls used to generate plot
        The available method names usually start with 'plot_'
      xlabel = x-axis label, as string
      ylabel = y-axis label, as string
      title = plot title, as string
      fmts = list of format specifier strings
    To be created by methods:
      series = sequence of PlotSeries instances
      fig = matplotlib Figure for the generated plot
      ax = matplotlib Axes for the plot
      info = dictionary of miscellaneous data"""
  __slots__=['figsize','filename','prepfunctions','plotfunctions','series','xlabel','ylabel','title','fmts','fig','ax','info']
  
  def execute_commandseq(self,attrname):
    """Execute the command sequence
    Arguments:
      attrname = name of attribute containing the command sequence"""
    for cmd in getattr(self,attrname,[]):
      #Function name and arguments
      funcname, kwargs = cmd
      #Call it
      try:
        getattr(self,funcname)(**kwargs)
      except Exception as einst:
        print("Excption occured for command: %s"%str(cmd), file=sys.stderr)
        raise einst
      
  def make_plot(self,*datafiles):
    """Create the plot.
    Arguments:
      *datafiles = data files, passed through to load_data"""
    #Load the data we need to generate the plot
    self.load_data(*datafiles)
    
    #Initialize the figure at the size requested
    self.fig = plt.figure(figsize=self.figsize)
    
    #Get the axes
    self.ax=self.fig.gca()
    
    #Call the preparation functions
    self.execute_commandseq('prepfunctions')
    
    #Add the available series to the axes
    self.plot_basic_series()
    
    #Call the requested plot functions
    self.execute_commandseq('plotfunctions')
    
    #Save the figure
    fpath=osp.join(self.outdir(),self.filename)
    os.makedirs(self.outdir(),exist_ok=True)
    self.fig.savefig(fpath)
    
    #Done
    return

  def plot_basic_series(self):
    """A simple plot."""
    for i,sr in enumerate(self.series):
        o=sr.add_to_axes(self.ax,self.fmts[i])
    if getattr(self,'title',None) is not None:
      o=self.ax.set_title(self.title)
    # if getattr(self,'legendloc',None) is not None:
    #   o=self.ax.legend(loc=self.legendloc)
    if getattr(self,'xlabel',None) is not None:
      o=self.ax.set_xlabel(self.xlabel)
    if getattr(self,'ylabel',None) is not None:
      o=self.ax.set_ylabel(self.ylabel)
    return   

  def plot_axmethod(self,method,kwargs=None):
    """Call a method of the axes.
    Arguments:
      method = name of Axes method to call, as string
      kwargs = arguments dictionary for the method"""
    f=getattr(self.ax,method)
    if kwargs is None:
      kwargs = {}
    f(**kwargs)
    return

  def plot_hline(self,locspec,kwargs=None):
    """Add a horizontal line with a value from info
    Arguments:
      locspec = sequence of keys in the info dictionary to locate the y-value
      kwargs = keyword arguments for ax.axhline"""
    yval=useful.nested_location(self.info,locspec)
    if kwargs is None:
      kwargs = {}
    self.ax.axhline(yval,**kwargs)
    return

  def plot_vline(self,locspec,kwargs=None):
    """Add a vertical line with a value from info
    Arguments:
      locspec = sequence of keys in the info dictionary to locate the x-value
      kwargs = keyword arguments for ax.axvline"""
    xval=useful.nested_location(self.info,locspec)
    if kwargs is None:
      kwargs = {}
    self.ax.axvline(xval,**kwargs)
    return


class ModelPlotFigure(PlotFigure):
  """Data for a single model plot
  Attributes:
    To be read in from yaml file:
      plotname = plot name in outdata file holding data series
    To be assigned after instantiation:
      modelname = name of model
    To be created by methods:
      (none)"""
  __slots__=['plotname','modelname']
  def outdir(self):
    return osp.join(postprocfolder,self.basename,self.modelname)

  def load_data(self,pklfile,infofile):
    """Load the data for the plot.
    Arguments:
      pklfile = path to the solver_general.OutData pickle file
      infofile = path to the info.yaml file"""
    #Load the data series
    outdata=solver_general.OutData.from_pickle(pklfile)
    self.series=outdata.plots[self.plotname]
    
    #Load the info
    self.info=useful.readyaml(infofile)
    
    return

  def prep_porelimits(self):
    self.info['meshparams']['poretop_z']=self.info['meshparams']['H']+self.info['meshparams']['tm']

class CollectionPlotFigure(PlotFigure):
  """Data for a single collection plot
  Attributes:
    To be read in from yaml file:
      calcfunctions = sequence of calculation functions to be called before generating plot
      seriescols = sequence of series definitions (xcol, ycol, label),
        where the columns specify the DataFrame columns containing values for the series
        The label is optional.
    To be created by methods:
      df = the DataFrame"""
  __slots__=['calcfunctions','seriescols','df']
  def outdir(self):
    return osp.join(postprocfolder,self.basename)

  def load_data(self,dfpath):
    """Load the data for the plot.
    Arguments:
      dfpath = path to the Pandas DataFrame to load"""
    #Load the DataFrame
    self.df=pd.read_pickle(dfpath)
    
    #Initialize empty info
    self.info={}
    
    #Do the requested calculations to add new columns
    self.execute_commandseq('calcfunctions')
    
    #Add the requested columns in as series
    self.series=[]
    for sdef in self.seriescols:
      sdef_dict={'xvals':self.df[sdef[0]],'yvals':self.df[sdef[1]]}
      if len(sdef)>2:
        sdef_dict['label']=sdef[2]
      else:
        sdef_dict['label']=''
      self.series.append(PlotSeries(**sdef_dict))
    
    return

  def calc_Dratio(self):
    def calc_ratio(row):
      return row['Deff']/row['D_bulk']
    self.df['ratio_D']=self.df.apply(calc_ratio,axis=1)
    return

  def prep_series_equality(self):
    pdser=self.df['free_volume_frac']
    vals=[pdser.min(),pdser.max()]
    ser=PlotSeries(xvals=vals,yvals=vals,label="1:1")
    self.series.append(ser)
    return