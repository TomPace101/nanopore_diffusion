#Apply post-processing defintions

#Standard library
import os.path as osp

#Site packages

#Local
import folderstructure
import solver_general
import plotdata
import collect_results
import useful

#Constants

class PostProcParameters(useful.ParameterSet):
  """Definition of post-processing requests
  Attributes:
    modelparamsfile = name of yaml file containing the models to post-process (include .yaml extension)
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
  __slots__=['modelparamsfile','modelnames','do_collection','collection_exclusions','model_plots','collection_plots']

def get_df_fname(basename):
  return osp.join(folderstructure.postprocfolder,basename,collect_results.collected_df_fname)

loaded_modelfiles=[]
loaded_models={}
def do_postproc(postprocparams):
  #Load the model parameters file if not already loaded
  if not postprocparams.modelparamsfile in loaded_modelfiles:
    modelparams_fpath=osp.join(folderstructure.params_model_folder,postprocparams.modelparamsfile)
    modelparams_gen=solver_general.ModelParameters.all_from_yaml(modelparams_fpath)
    for mp in modelparams_gen:
      loaded_models[mp.modelname]=mp
    loaded_modelfiles.append(postprocparams.modelparamsfile)
  #Get the indicated modelnames
  if getattr(postprocparams,'modelnames',None) is not None:
    modellist=[v for k,v in loaded_models.items() if k in postprocparams.modelnames]
  else:
    modellist=[m for m in loaded_models.values()]
  #Generate model plots, if any
  if getattr(postprocparams,'model_plots',None) is not None:
    #Each requested model plot
    for mplotdict in postprocparams.model_plots:
      for modelparams in modellist:
        mplot = plotdata.ModelPlotFigure(**mplotdict)
        mplot.basename=postprocparams.basename
        mplot.modelname=modelparams.modelname
        mplot.make_plot()
  #Do the collection, if requested
  if getattr(postprocparams,'do_collection',False):
    inputfiles=collect_results.list_inputfiles(postprocparams.basename,modellist)
    outfpath=get_df_fname(postprocparams.basename)
    exclusions=getattr(postprocparams,'collection_exclusions',[])
    collect_results.do_collection(inputfiles,outfpath,exclusions)
    #Generate each requested plot from the collected data
    if getattr(postprocparams,'collection_plots',None) is not None:
      for cplotdict in postprocparams.collection_plots:
        cplot=plotdata.CollectionPlotFigure(**cplotdict)
        cplot.basename=postprocparams.basename
        cplot.make_plot()

#Support command-line arguments
if __name__ == '__main__':
  program_description='Run post-processing commands'
  input_file_description="Path to the file containing PostProcParameters definitions"
  
  useful.run_cmd_line(program_description,input_file_description,PostProcParameters,do_postproc)
