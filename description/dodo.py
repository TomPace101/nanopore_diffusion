#Doit file for problem description

import os
import os.path as osp

#Constants
svgdir='fig_svg'
pdfdir='fig_pdf'
latexinputs=['resultcalcs.tex','eqn_uh_fick.tex','geometry.tex','eqn_uh_smol.tex']

#Get the names of all the image files
def list_svg_stems(svgdir):
  "get a list of stems for the svg files"
  allfiles=os.listdir(svgdir)
  allsplits=[osp.splitext(x) for x in allfiles]
  return [stem for stem,ext in allsplits if ext.lower()=='.svg']

def stem_to_file(stems,folder,ext):
  "get file names with a folder and extension from a list of stem names"
  return [osp.join(folder,'%s.%s'%(x,ext)) for x in stems]

figstems=list_svg_stems(svgdir)
fig_inputs=stem_to_file(figstems,svgdir,'svg')
fig_outputs=stem_to_file(figstems,pdfdir,'pdf')

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
  