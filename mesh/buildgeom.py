#Process the jinja2 template into one (or more) gmsh .geo file(s)
#Usage:
#python buildgeom.py geomdef paramdef
#for more details, see
#python buildgeom.py -h

import argparse
import os.path as osp

from jinja2 import Template
import yaml

#Process command-line arguments
parser = argparse.ArgumentParser(description='Create gmsh .geo file')
parser.add_argument('geomdef', help="geometry definition yaml file, which provides the template data needed by the template")
parser.add_argument('paramdef', help="parameter defnition yaml file, which assigns values to select parameters in the template")
cmdline=parser.parse_args()
assert osp.isfile(cmdline.geomdef), "Geometry definition file does not exist: %s"%cmdline.geomdef
assert osp.isfile(cmdline.paramdef), "Parameter definition file does not exist: %s"%cmdline.paramdef

#Read in the two yaml files
with open(cmdline.geomdef,'r') as fp:
  dat=fp.read()
  geomdef=yaml.load(dat)
with open(cmdline.paramdef,'r') as fp:
  dat=fp.read()
  paramdef=yaml.load(dat)

#Namepsace the geometry definition for convenience
geom=argparse.Namespace(**geomdef)

#Put needed parameters into template input
t_input={}
t_input.update(paramdef)
t_input['mshfile']=geom.mshfile
t_input['ptstrs']=[str(x) for x in geom.ptlist]

#From mapping of surfaces to points, generate:
# - mapping of loops to line and circle names
# - mapping of line and circle names to list of points
loops={}
lines={}
circles={}
lnum=1
cnum=1
for surfnum, pttup in geom.geomtable.items():
  loops[surfnum]=[]
  startpt=pttup[0]
  indx=1
  while indx < len(pttup):
    if pttup[indx]=='center':
      indx += 2
      cpoint=pttup[indx-1]
      endpt=pttup[indx]
      ctup=(startpt,cpoint,endpt)
      rctup=(endpt,cpoint,startpt) #Reverse order
      #Already exists?
      found=False
      for cn, cpts in circles.items():
        if cpts==ctup:
          found=True
          loops[surfnum].append(cn)
          break
        elif cpts==rctup:
          found=True
          loops[surfnum].append('-'+cn)
          break
      if not found:
        cname='c%d_%d_%d'%ctup
        circles[cname]=ctup
        loops[surfnum].append(cname)
        cnum += 1
    else:
      endpt=pttup[indx]
      ltup=(startpt,endpt)
      rltup=(endpt,startpt) #Reverse order
      #Already exists?
      found=False
      for ln, lpts in lines.items():
        if lpts==ltup:
          found=True
          loops[surfnum].append(ln)
          break
        elif lpts==rltup:
          found=True
          loops[surfnum].append('-'+ln)
          break
      if not found:
        lname='l%d_%d'%ltup
        lines[lname]=ltup
        loops[surfnum].append(lname)
        lnum += 1
    #Next point
    startpt=pttup[indx]
    indx += 1
#Provide mappings to template
linemap=dict([(n,', '.join(['p%d'%p for p in pts])) for n,pts in lines.items()])
t_input['lines']=linemap
circmap=dict([(n,', '.join(['p%d'%p for p in pts])) for n, pts in circles.items()])
t_input['circles']=circmap
loopmap=dict([(n,', '.join([x for x in ents])) for n,ents in loops.items()])
t_input['loops']=loopmap
#Apply surface types
surftypes=dict([(x,'Ruled' if x in geom.nonplanar else 'Plane') for x in geom.geomtable.keys()])
t_input['surftypes']=surftypes
#Apply reversal to selected surfaces
surfnums=[-x if x in geom.revsurfs else x for x in geom.geomtable.keys()]
t_input['looplist']=', '.join([str(x) for x in surfnums])

#Read template
with open(geom.tmplfile,'r') as fp:
    tmpldat=fp.read()

#Render template
tmplobj=Template(tmpldat, trim_blocks=True)
outdat = tmplobj.render(t_input)

#Output result
with open(geom.outfile,'w') as fp:
    fp.write(outdat)
