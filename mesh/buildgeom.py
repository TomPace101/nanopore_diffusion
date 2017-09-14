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
t_input['ptstrs']=dict([(str(x),y) for x,y in geom.ptdict.items()])

#From mapping of surfaces to points, generate:
# - mapping of loops to line and circle names
# - mapping of line and circle names to list of points
def add_entity(tup, tdict, looplist, nameprefix):
  """Create lines/circles needed for line loops, unless they already exist.
  Inputs:
    tup = line/circle tuple of point numbers (as integers) (lines have 2 points, circles have 3)
    tdict = dictionary storing the line/circle tuples
    looplist = list of lines/circles in this line loop
    nameprefix = 'c' for circle, 'l' for line
  No return value.
  Side effects:
    tdict and looplist are modified in place, to include the generated (or found) lines/circles"""
  #Reversed tuple
  rtup = tuple(reversed(tup))
  found=False
  #Search through all lines/circles already created, for one with same set of points, in forward or reverse order
  for n, pts in tdict.items():
    if pts==tup:
      found=True
      looplist.append(n)
      break
    elif pts==rtup:
      found=True
      looplist.append('-'+n)
      break
  #Create new line/circle if not found
  if not found:
    nametmpl='_'.join('%d' for i in range(len(tup)))
    name=nameprefix+nametmpl%tup
    tdict[name]=tup
    looplist.append(name)
  return

#Create dictionaries of lines, circles, and line loops
loops={}
lines={}
circles={}
for surfnum, pttup in geom.geomtable.items():
  loops[surfnum]=[]
  startpt=pttup[0]
  indx=1
  while indx < len(pttup):
    if pttup[indx]=='center':
      indx += 2
      ctup=(startpt,pttup[indx-1],pttup[indx])
      add_entity(ctup,circles,loops[surfnum],'c')
    else:
      ltup=(startpt,pttup[indx])
      add_entity(ltup,lines,loops[surfnum],'c')
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

#Apply reversal to selected surfaces for surface loops
surfnums=[-x if x in geom.revsurfs else x for x in geom.geomtable.keys()]
t_input['looplist']=', '.join([str(x) for x in surfnums])

#Read template
with open(geom.tmplfile,'r') as fp:
    tmpldat=fp.read()

#Render template
tmplobj=Template(tmpldat, trim_blocks=True)
outdat = tmplobj.render(t_input)

#Output result
with open(paramdef['outfile'],'w') as fp:
    fp.write(outdat)
