#Process the jinja2 template into one (or more) gmsh .geo file(s)
#Usage:
#python buildgeom.py

from jinja2 import Template

## TODO: read parameters in from command line, or yaml file

#Just dummy parameters for now
params={'mcar':0.01, 'Lx':1.5, 'Ly':1.0, 'R':0.75, 'H':1.0, 'tm':2.0}

#Hardcoded file paths
## TODO: de-hardcode
tmplfile='body-centered.geo.jinja2'
outfile='body-centered.geo'

#Hardcoded geometry information
## TODO: de-hardcode (should be associated with jinja2 file)
#Point list
ptlist = [111,      311, 331, 131,
          112, 212, 312, 332, 132, 122,
          113, 213, 313, 333, 133, 123,
          114,      314, 334, 134]
#Point list, as strings, is needed by template
params['ptstrs']=[str(x) for x in ptlist]
#Mapping of surfaces to points
#In the following table of values, each tuple is a line loop that will be turned into a surface.
#Circles can be included by including 'center' followed by the center point.
#The preceeding and next points are the circle start and end points, respectively.
geomtable = {1:(111, 311, 331, 131, 111),
             2:(212, 312, 332, 132, 122, 'center', 112, 212),
             3:(213, 313, 333, 133, 123, 'center', 113, 213),
             4:(114, 314, 334, 134, 114),
             5:(111, 311, 312, 212, 213, 313, 314, 114, 111),
             6:(111, 131, 132, 122, 123, 133, 134, 114, 111),
             7:(331, 131, 132, 332, 331),
             8:(311, 331, 332, 312, 311),
             9:(333, 133, 134, 334, 333),
             10:(313, 333, 334, 314, 313),
             11:(212, 'center', 112, 122, 123, 'center', 113, 213, 212)}

#From mapping of surfaces to points, generate:
# - mapping of loops to line and circle names
# - mapping of line and circle names to list of points
loops={}
lines={}
circles={}
lnum=1
cnum=1
for surfnum, pttup in geomtable.items():
  loops[surfnum]=[]
  startpt=pttup[0]
  indx=1
  while indx < len(pttup):
    if pttup[indx]=='center':
      indx += 2
      cpoint=pttup[indx-1]
      endpt=pttup[indx]
      circles[cnum]=(startpt,cpoint,endpt)
      loops[surfnum].append('c%d'%cnum)
      cnum += 1
    else:
      endpt=pttup[indx]
      lines[lnum]=(startpt,endpt)
      loops[surfnum].append('l%d'%lnum)
      lnum += 1
    #Next point
    startpt=pttup[indx]
    indx += 1
#Provide mappings to template
linemap=dict([('l%d'%n,', '.join(['p%d'%p for p in pts])) for n,pts in lines.items()])
params['lines']=linemap
circmap=dict([('c%d'%n,', '.join(['p%d'%p for p in pts])) for n, pts in circles.items()])
params['circles']=circmap


#Read template
with open(tmplfile,'r') as fp:
    tmpldat=fp.read()

#Render template
tmplobj=Template(tmpldat, trim_blocks=True)
outdat = tmplobj.render(params)

#Output result
with open(outfile,'w') as fp:
    fp.write(outdat)
