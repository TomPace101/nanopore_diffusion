#Doit file for problem description

import os
import os.path as osp

#Constants
svgdir='fig_svg'
pdfdir='fig_pdf'

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
  infile='description.md'
  outfile='description.pdf'
  return {'actions': ['pandoc --number-sections -o %s %s'%(outfile,infile)],
          'file_dep': [infile]+fig_outputs,
          'targets': [outfile]}

def task_pdf_svg():
  filepairs=zip(fig_inputs,fig_outputs)
  for infile,outfile in filepairs:
    tdef={'actions':['inkscape --without-gui --export-area-drawing --export-pdf=%s %s'%(outfile,infile)],
          'file_dep':[infile],
          'targets':[outfile],
          'name':outfile}
    yield tdef
  