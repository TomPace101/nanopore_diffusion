#Doit file for problem description

import os
import os.path as osp

#Constants
svgdir='fig_svg'
pdfdir='fig_pdf'
partsdir='parts'

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

latexinputs=list_files(partsdir,'tex')

def task_gen_report():
  stemname='description'
  infile=stemname+'.tex'
  outfile=stemname+'.pdf'
  cleancmd = '; '.join(['rm -f %s.%s'%(stemname,ext) for ext in ['aux','log','out','pdf']])
  cmdstr='pdflatex -interaction=nonstopmode -halt-on-error %s'%(infile)
  return {'actions': [cleancmd, cmdstr, cmdstr], #clean up first, then run pdflatex twice to get figure references correct.
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
  