"""Support for building FEniCS weak forms"""

#Standard library
from collections import OrderedDict

class EquationTerm(object):
  """To keep track of additional information about each UFL form added up into the equation to solve.
  
  Attributes:
  
    - name = a string uniquely identifying the equation term
    - ulf = the UFL form object for the term"""

  def __init__(self,name,ufl,**kwargs):
    self.name=name
    self.ufl=ufl
    for k,v in kwargs.items():
      setattr(self,k,v)

class EquationTermDict(OrderedDict):
  """An ordered dictionary of equation terms.
  
  Attributes:
  
    - termclass = the class of which all equationterms are instances"""

  def __init__(self,termclass=EquationTerm,*args,**kwargs):
    #Initialization from base class
    super(EquationTermDict, self).__init__(*args,**kwargs)
    #Store the term class
    self.termclass=termclass

  def add(self,*args,**kwargs):
    """Create a new term and add it to this dictionary.
    
    Arguments are the same as those to create a new EquationTerm"""
    term=self.termclass(*args,**kwargs)
    self[term.name]=term

  def selectterms(self,**kwargs):
    """Return an EquationTermDict of terms with the specified properties"""
    out=EquationTermDict(self.termclass)
    for name,term in self.items():
      include=True
      for k,v in kwargs.items():
        if getattr(term,k)!=v:
          include=False
          break
      if include:
        out[name]=term
    return out

  def sumterms(self,zeroval=None,**kwargs):
    """Return the sum of terms with the specified properties, as a UFL form
    
    Arguments:
    
      - zeroval = optional, UFL form to be returned if no terms match the specified properties
      - **kwargs = keyword arguments for selectterms"""
    terms=self.selectterms(**kwargs)
    if len(terms)>0:
      return sum([t.ufl for t in terms.values()])
    elif zeroval is not None:
      return zeroval
    else:
      #zeroval was not provided
      raise RuntimeError("Unable to return a proper sum of zero terms; Provide keyword argument zeroval to resolve.")
