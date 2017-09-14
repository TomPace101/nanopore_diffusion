
def task_gen_report():
  infile='description.md'
  outfile='description.pdf'
  return {'actions': ['pandoc --number-sections -o %s %s'%(outfile,infile)],
          'file_dep': [infile],
          'targets': [outfile]}