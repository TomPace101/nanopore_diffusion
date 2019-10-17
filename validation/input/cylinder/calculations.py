"""Post-processing calculations"""

import numpy as np

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
  
    - freevol = optional, column name for the radius
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


#List of functions to be bound as methods
request_methods=[calc_porosity_column_byradius]
