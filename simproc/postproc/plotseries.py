"""Defining the data series that go on plots"""


class PlotSeries:
  """Data for a single series on a plot

  Attributes:

    - xvals = array of x-values
    - yvals = array of y-values
    - label = legend label
    - metadata = other parameters needed to identify the data series
  """
  def __init__(self,xvals,yvals,label=None,metadata=None):
    self.xvals=xvals
    self.yvals=yvals
    self.label=label
    self.metadata=metadata
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

##Should this be a classmethod of PlotSeries instead of a class?
#argue no: not loadable from yaml then
#argument yes: don't need to load from yaml
class PlotDataFrameSeries:
  """Define a PlotSeries from a DataFrame"""
  def __new__(cls,df,xcol,ycol,query=None,label=None,metadata=None):
    #Select particular rows if requested
    if query is None or len(query)==0:
      qdf=df
    else:
      qdf=df.query(query)
    #Get the data for the series from the dataframe
    xvals=qdf[xcol]
    yvals=qdf[ycol]
    #Instantiate!
    return PlotSeries(xvals=vals,yvals=yvals,label=label,metadata=metadata)
