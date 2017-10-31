#Doit file for problem description

import os
import os.path as osp
import shutil

#Constants
svgdir='fig_svg'
pdfdir='fig_pdf'
partsdir='parts'
outdir='output'
post_clean_exts=['aux','log','out','bbl','bcf','blg','run.xml','toc']
pre_clean_exts=post_clean_exts+['pdf']

def list_stems(folder,ext):
  "get a list of stems for the files in the directory with the given extension"
  allfiles=os.listdir(folder)
  allsplits=[osp.splitext(x) for x in allfiles]
  return [stem for stem,xt in allsplits if xt.lower()[1:]==ext]

def stem_to_file(stems,folder,ext):
  "get file names with a folder and extension from a list of stem names"
  return [osp.join(folder,'%s.%s'%(x,ext)) for x in stems]

def list_files(sdir,ext):
  "get a list of the files in a directory with the specified extension"
  stems=list_stems(sdir,ext)
  return stem_to_file(stems,sdir,ext)

figstems=list_stems(svgdir,'svg')
fig_inputs=stem_to_file(figstems,svgdir,'svg')
fig_outputs=stem_to_file(figstems,pdfdir,'pdf')

latexinputs=list_files(partsdir,'tex')+list_files(partsdir,'bib')

def do_post_clean(all_list):
  "Move the multitude of latex-generated output files to a subfolder"
  for fname in all_list:
    if osp.isfile(fname):
      newfile=osp.join(outdir,fname)
      if osp.isfile(newfile):
        os.remove(newfile)
      shutil.move(fname,osp.join(outdir,fname))

def task_gen_report():
  stemname='description'
  infile=stemname+'.tex'
  outfile=stemname+'.pdf'
  pre_clean_cmd = '; '.join(['rm -f %s.%s'%(stemname,ext) for ext in pre_clean_exts])
  post_clean_files=['%s.%s'%(stemname,ext) for ext in post_clean_exts]
  cmdstr='pdflatex -interaction=nonstopmode -halt-on-error %s'%(infile)
  bibcmd='biber '+stemname
  return {'actions': [pre_clean_cmd, cmdstr, bibcmd, cmdstr, (do_post_clean,(post_clean_files,))], #clean up first, then run pdflatex twice to get figure references correct, with bibtex in between.
          'file_dep': [infile]+fig_outputs+latexinputs,
          'targets': [outfile]}

def task_pdf_svg():
  filepairs=zip(fig_inputs,fig_outputs)
  for infile,outfile in filepairs:
    tdef={'actions':['inkscape --without-gui --export-area-drawing --export-pdf=%s %s'%(outfile,infile)],
          'file_dep':[infile],
          'targets':[outfile],
          'name':outfile}
    yield tdef
  