#Process the jinja2 template into one (or more) gmsh .geo file(s)
#Usage:
#python buildgeom.py

from jinja2 import Template

## TODO: read parameters in from command line, or yaml file

#Just dummy parameters for now
params={'mcar':0.01, 'Lx':1.5, 'Ly':1.0, 'R':0.75, 'H':1.0, 'tm':2.0}

#Hardcoded file paths
tmplfile='body-centered.geo.jinja2'
outfile='body-centered.geo'


#Read template
with open(tmplfile,'r') as fp:
    tmpldat=fp.read()

#Render template
tmplobj=Template(tmpldat)
outdat = tmplobj.render(params)

#Output result
with open(outfile,'w') as fp:
    fp.write(outdat)
