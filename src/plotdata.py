
#Standard library
import os
import os.path as osp

#Site packages
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Local
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
      plotfunction = name of method used to generate plot
      xlabel = x-axis label, as string
      ylabel = y-axis label, as string
      title = plot title, as string
      fmts = list of format specifier strings
      legendloc = legend location (loc keyword argument for Axes.legend()), None or omitted for no legend
    To be created by methods:
      fig = matplotlib Figure for the generated plot
      ax = matplotlib Axes for the plot"""
  __slots__=['figsize','filename','plotfunction','xlabel','ylabel','title','fmts','legendloc','fig','ax']
  def make_plot(self,*datafiles):
    """Create the plot.
    Arguments:
      *datafiles = data files, passe through to load_data"""
    #Load the data we need to generate the plot
    self.load_data(*datafiles)
    
    #Initialize the figure at the size requested
    self.fig = plt.Figure(figsize=self.figsize)
    
    #Get the axes
    self.ax=fig.gca()
    
    #Call the requested plot function
    plotfunc=getattr(self,self.plotfunction)
    plotfunc()
    
    #Save the figure
    os.makedirs(osp.split(self.filename)[0],exist_ok=True)
    self.fig.savefig(self.filename)
    
    #Done
    return

  def plot_basic_series(self):
    """A simple plot."""
    for i,sr in enumerate(self.series):
        o=sr.add_to_axes(self.ax,self.fmts[i])
    if getattr(self,'title',None) is not None:
      o=self.ax.set_title(self.title)
    if getattr(self,'legendloc',None) is not None:
      o=self.ax.legend(loc=self.legendloc)
    if getattr(self,'xlabel',None) is not None:
      o=self.ax.set_xlabel(self.xlabel)
    if getattr(self,'ylabel',None) is not None:
      o=self.ax.set_ylabel(self.ylabel)
    return   


class ModelPlotFigure(PlotFigure):
  """Data for a single model plot
  Attributes:
    To be read in from yaml file:
      plotname = plot name in outdata file holding data series
    To be created by methods:
      series = sequence of PlotSeries instances"""
  __slots__=['plotname',]
  def load_data(self,pklfile,infofile):
    """Load the data for the plot.
    Arguments:
      pklfile = path to the solver_general.OutData pickle file
      infofile = path to the info.yaml file"""
    #Load the data series
    oudata=solver_general.OutData.from_pickle(pklfile)
    self.series=outdata.plots[self.plotname]
    
    #Load the info
    self.info=useful.readyaml(infofile)
    
    return

  def plot_invert_xaxis(self):
    "just invert the x-axis"
    self.plot_basic_series()
    self.ax.invert_xaxis()
    return

  def plot_radial_potential(self):
    pore_radius = self.info['meshparams']['R']
    applied_potential=self.info['conditions']['potential']['bcdict'][11]
    self.plot_basic_series()
    o=self.ax.axvline(pore_radius,label='Pore Boundary',color='k',linestyle='--')
    o=self.ax.axhline(applied_potential,label='Potential at Pore Boundary',color='k',linestyle=":")
    

class CollectionPlotFigure(PlotFigure):
  """Data for a single collection plot
  Attributes:
    To be read in from yaml file:
      calcfuncs = sequence of calculation functions to be called before generating plot
      xcol = sequence of column names in DataFrame to use for x-values
      ycol = sequence of column names in DataFrame to use for y-values
    To be created by methods:
      df = the DataFrame"""
  __slots__=['calcfuncs','xcol','ycol','df']
  def load_data(self,dfpath):
    """Load the data for the plot.
    Arguments:
      dfpath = path to the Pandas DataFrame to load"""
    #Load the DataFrame
    self.df=pd.read_pickle(dfpath)
    
    #Do the requested calculations to add new columns
    if hasattr(self,'calcfuncs') and self.calcfuncs is not None:
      for cf in self.calcfuncs:
        f=getattr(self,cf)
        f()
    
    return

    def calc_Dratio(self):
      def calc_ratio(row):
        return row['Deff']/row['D_bulk']
      self.df['ratio_D']=self.df.apply(calc_ratio,axis=1)
      return