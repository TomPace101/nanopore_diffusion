#Apply post-processing defintions

#Standard library

#Site packages

#Local
from folderstructure import *
import useful

#Constants

class PostProcParameters(useful.ParameterSet):
  """Definition of post-processing requests
  Attributes:
    modelnames = sequence of model names that this post-processing definition applies to
      Omit or set to None to apply to all model names for the basename.
      An empty list will result in no post-processing being done.
    do_collection = boolean, True to collect info.yaml files into a Pandas DataFrame
    collection_exclusions = list of keys in info.yaml to be excluded from DataFrame
    model_plots = sequence of plots to generate for each model
      Each plot defines a plotdata.PlotFigure instance.
    collection_plots = sequence of plots to generate from the collected DataFrame
      Each plot defines a plotdata.PlotFigure instance.
      If collection_schema is omitted or None, this should be as well, or be an empty sequence."""
  __slots__=['modelnames','do_collection','collection_exclusions','model_plots','collection_plots']