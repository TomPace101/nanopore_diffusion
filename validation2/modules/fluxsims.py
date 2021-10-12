"""Pre/Post-processing calculations for the direct Smoluchowski simulations"""

import fenics as fem

def effective_D(self,outattr,fluxattr,areaattr,start_conc_attr,end_conc_attr,delta_s_attr):
  """Calculate effective diffusion constant, using concentration averages rather than point values

  Arguments:

    - outattr = attribute path for storage of result
    - fluxattr = attribute path to previously calculated total

        This requires a previous call to fluxintegral (or similar).

    - areaattr = attribute path to previously calculated area in results dictionary

        This requires a previous call to facet_area (or similar).

    - start_conc_attr = attribute path to concentration value at starting boundary
    - end_conc_attr = attribute path to concentration value at ending boundary
    - delta_s_attr = attribute path to value of Delta s for effective D calculation

  No return value.
  No output files."""
  #Get the values from attributes
  totflux=self.get_nested(fluxattr)
  area=self.get_nested(areaattr)
  startconc=self.get_nested(start_conc_attr)
  endconc=self.get_nested(end_conc_attr)
  delta_s=self.get_nested(delta_s_attr)
  #Calculate the change in concentration between the two points
  delta_c=endconc-startconc
  #Calculate diffusion constant
  Deff=float(totflux/area*delta_s/delta_c)
  #Store result
  self.set_nested(outattr,Deff)
  return

#List of functions to be bound as methods
request_methods=[effective_D]
