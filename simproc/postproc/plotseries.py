
class PlotSeries:
  """Data for a single series on a plot

  Attributes:

    - xvals = array of x-values
    - yvals = array of y-values
    - label = legend label
    - metadata = other parameters needed to identify the data series
  """
  _expected_attrs=['xvals','yvals','label','metadata']
  def __init__(self,**kwargs):
    for k,v in kwargs.items():
      setattr(self,k,v)
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
      label=self.label
    else:
      label=newlabel
    return ax.plot(self.xvals,self.yvals,fmt,label=label,**kwd)
