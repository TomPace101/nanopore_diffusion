"""Module holding methods that are used by more than one simulator subclass, but not enough to be in simrequest itself."""

import fenics as fem

def calcflux(self, solnattr='soln', idx=None, attrpath='flux', Dattr='Dbar_proj'):
  """Flux as vector field (new attribute)"""
  soln=self.get_nested(solnattr)
  Dvalue=self.get_nested(Dattr)
  funcname="flux_"+solnattr
  if idx is not None:
    soln=soln[idx]
    Dvalue=Dvalue[idx]
    funcname += "_%d"%idx
  expr=-Dvalue*fem.grad(soln)
  fluxres=fem.project(expr,self.V_vec,solver_type="cg",preconditioner_type="amg") #Solver and preconditioner selected to avoid UMFPACK "out of memory" error (even when there's plenty of memory)
  fluxres.rename(funcname,"calculated flux")
  self.set_nested(attrpath,fluxres)
  return

def fluxintegral(self,pfacet,attrpath,internal=False,fluxsign=None,normalvar=None,fluxattr='flux',idx=None):
  """Flux integral over specified facet

  Arguments:
    pfacet = physical facet number for flux measurement
    attrpath = attribute path for storage of the integral
    internal = boolean, default False, True to use internal boundary, False for external
    fluxsign = '+' or '-' to specify which direction normal to the facet for flux calculation
      Required only if internal==True
    fluxattr = name of attribute containing the flux field for integration
  No return value.
  No output files."""
  n=fem.FacetNormal(self.meshinfo.mesh)
  if internal:
    integral_type='interior_facet'
    assert fluxsign=='+' or fluxsign=='-', "Invalid fluxsign: %s"%str(fluxsign)
    this_n=n(fluxsign)
  else:
    integral_type='exterior_facet'
    this_n=n
  this_ds=fem.Measure(integral_type,domain=self.meshinfo.mesh,subdomain_data=self.meshinfo.facets)
  totflux=fem.assemble(fem.dot(getattr(self,fluxattr),this_n)*this_ds(pfacet))
  self.set_nested(attrpath,totflux)
  return
