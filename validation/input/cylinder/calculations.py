"""Post-processing calculations"""

import numpy as np

from simproc.postproc.plotseries import PlotSeries

def calc_porosity_column(self, freevol, totvol, newcol="porosity", dfpath="df"):
  """Compute the porosity column of the dataframe

  Arguments:
  
    - freevol = column name for the free volume
    - totvol = column name for the total volume
    - newcol = optional, name of new column storing results, as string
    - dfpath = optional, attribute path to the dataframe

  New column added to the dataframe.
  No return value.
  No output files."""
  df=self.get_nested(dfpath)
  def calcphi(row):
    return row[freevol]/row[totvol]
  df[newcol]=df.apply(calcphi,axis=1)
  return

def calc_porosity_column_byradius(self, radcol="r", newcol="porosity", dfpath="df"):
  """Compute the porosity column of the dataframe

  Arguments:
  
    - freevol = optional, column name for the radius, as string
    - newcol = optional, name of new column storing results, as string
    - dfpath = optional, attribute path to the dataframe

  New column added to the dataframe.
  No return value.
  No output files."""
  df=self.get_nested(dfpath)
  def calcphi(row):
    r=row[radcol]
    return 1-np.pi*r**2
  df[newcol]=df.apply(calcphi,axis=1)
  return

def calc_cylinder_hashin_shtrikman(self, phicol="porosity", newcol="upper_bound", dfpath="df"):
  """Compute the Hashin-Shtrikman upper bound for the cylinder problem

  Arguments:

    - phicol = optional, column name for the porosity, as string
    - newcol = optional, column name for storing results, as string
    - dfpath = optional, attribute path to the dataframe

  New column added to the dataframe.
  No return value.
  No output files."""
  df=self.get_nested(dfpath)
  def calcval(row):
    phi=row[phicol]
    return 2*phi/(3-phi)
  df[newcol]=df.apply(calcval,axis=1)
  return

def series_hashin_shtrikman(self,start=0.2,stop=1.0,numpts=100,attrpath="hs_ub",label=None,metadata=None):
  """Compute a series of Hashin-Shtrikman upper bounds
  (As opposed to only computing them at simulation porosities.)

  Arguments:

    - start = optional, starting porosity value
    - stop = optional, final porosity value
    - numpts = optional, number of curve points
    - attrpath = optional, attribute path for storing results, as string
      - label = optional series label
      - metadata = optional series metadata

  New attribute added/overwritten.
  No return value.
  No output files."""

  phivals=np.linspace(start,stop,numpts)
  yvals=np.array([2*phi/(3-phi) for phi in phivals])
  series=PlotSeries(xvals=phivals,yvals=yvals,label=label,metadata=metadata)
  self.set_nested(attrpath,series)
  return


#List of functions to be bound as methods
request_methods=[calc_porosity_column_byradius, calc_cylinder_hashin_shtrikman, series_hashin_shtrikman]
