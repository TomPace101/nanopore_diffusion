
#Standard library
import os
import os.path as osp

#Site packages
from doit.tools import config_changed

#Local
from folderstructure import *
import useful
import collect_results
import plotdata
import postproc

#Constants
collected_df_fname='collected_results.pkl.gz'

def get_df_fname(basename):
  return osp.join(postprocfolder,basename,collected_df_fname)

def do_collection(basename,model_list,exclusions):
  infiles=collect_results.list_inputfiles(basename,model_list)
  outfpath=get_df_fname(basename)
  tdef = {'name': basename+":collection",
          'file_dep':infiles+[collect_results.__file__],
          'targets':[outfpath],
          'actions':[(collect_results.do_collection,(infiles,outfpath,exclusions))]}
  return tdef

def do_collectionplot(cplotparams,basename):
  dfpath=get_df_fname(basename)
  outfpath=osp.join(postprocfolder,basename,cplotparams.filename)
  codefile=osp.join(srcfolder,'plotdata.py')
  tdef = {'name':basename+":"+cplotparams.filename,
          'file_dep':[dfpath,codefile],
          'uptodate':[config_changed(cplotparams.to_dict())],
          'targets':[outfpath],
          'actions':[(cplotparams.make_plot,(dfpath,))]}
  return tdef

def do_modelplot(mplotparams,modelparams,basename):
  outfpath=osp.join(postprocfolder,basename,modelparams.modelname,mplotparams.filename)
  datadir=osp.join(solnfolder,basename,modelparams.modelname)
  codefile=osp.join(srcfolder,'plotdata.py')
  pklfile=osp.join(datadir,'outdata.pkl')
  infofile=osp.join(datadir,'info.yaml')
  tdef = {'name':modelparams.modelname+":"+mplotparams.filename,
          'file_dep':[pklfile,infofile,codefile],
          'uptodate':[config_changed(mplotparams.to_dict())],
          'targets':[outfpath],
          'actions':[(mplotparams.make_plot,(pklfile,infofile))]}
  return tdef

def postproc_task_generator(model_infiles,allmodels):
  "generator yielding postproc tasks"
  for fname in model_infiles:
    basename = osp.splitext(fname)[0]
    postproc_file = osp.join(params_postproc_folder,fname)
    postproc_list = [p for p in postproc.PostProcParameters.all_from_yaml(osp.join(params_postproc_folder,fname))]
    for postprocparams in postproc_list:
      #Get the indicated modelnames
      if getattr(postprocparams,'modelnames',None) is not None:
        modellist=[v for k,v in allmodels.items() if k in postprocparams.modelnames]
      else:
        modellist=[m for m in allmodels.values()]
      #Generate model plot tasks, if any
      if getattr(postprocparams,'model_plots',None) is not None:
        #Each requested model plot
        for mplotdict in postprocparams.model_plots:
          mplotparams = plotdata.ModelPlotFigure(**mplotdict)
          mplotparams.basename=postprocparams.basename
          for modelparams in modellist:
            mplotparams.modelname=modelparams.modelname
            yield do_modelplot(mplotparams,modelparams,basename)
      #Generate the collection task, if it exists
      if getattr(postprocparams,'do_collection',False):
        yield do_collection(basename,modellist,getattr(postprocparams,'collection_exclusions',[]))
        #Generate tasks for each plot from the collected data
        if getattr(postprocparams,'collection_plots',None) is not None:
          for cplotdict in postprocparams.collection_plots:
            cplotparams=plotdata.CollectionPlotFigure(**cplotdict)
            cplotparams.basename=postprocparams.basename
            yield do_collectionplot(cplotparams,basename)
  