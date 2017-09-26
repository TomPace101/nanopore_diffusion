#Process the jinja2 template into one (or more) gmsh .geo file(s)
#Usage:
#python buildgeom.py geomdef paramdef
#for more details, see
#python buildgeom.py -h

## TODO: validation of geometric inputs (different formulas for different geometries)

#Standard library
import argparse
import os.path as osp

#Site packages
from jinja2 import Environment, FileSystemLoader

#Local
import useful

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

def prepare_template_input(geom, paramdef):
  """Prepare the input dictionary for a template.
  Inputs: see write_one_geo, except that geom is a namespace rather than a dictionary, and outfile is not required
  Returns:
    t_input = the input dictionary for the template specified in geomdef"""

  #Put needed parameters into template input
  t_input={}
  t_input.update(paramdef)
  t_input['ptstrs']=dict([(str(x),y) for x,y in geom.ptdict.items()])

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
  
  return t_input

def write_one_geo(geomdef, paramdef, geofile):
  """Generate a single geo file based on a geometry defintion dictionary and parameter dictionary
  Inputs:
    geomdef = geometry definition dictionary, which must contain:
      tmplfile: geometry template file
      ptdict: dictionary of points and their corresponding mesh density parameter name
      geomtable: mapping of surfaces to points
      revsurfs: list of surfaces needing orientation reversal
      nonplanar: list of surfaces that are not planar surfaces
    paramdef = parameter defintion dictionary, which must contain:
      outfile: the .geo file to write
      mshfile: the .msh file for gmsh to create
      and all the other parameters needed by the geometry template file
    geofile = path to output .geo file, as string
  No return value. The .geo file is written."""

  #Namepsace the geometry definition for convenience
  geom=argparse.Namespace(**geomdef)

  #Get the input dictionary for the template
  t_input = prepare_template_input(geom, paramdef)

  #Load template
  env=Environment(loader=FileSystemLoader('.'), trim_blocks=True)
  tmpl=env.get_template(geom.tmplfile)

  #Render template
  outdat = tmpl.render(t_input)

  #Output result
  with open(geofile,'w') as fp:
    fp.write(outdat)
  return

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  parser = argparse.ArgumentParser(description='Create gmsh .geo file')
  parser.add_argument('geomdef', help="geometry definition yaml file, which provides the template data needed by the template")
  parser.add_argument('paramdef', help="parameter defnition yaml file, which assigns values to select parameters in the template")
  parser.add_argument('geofile', help="name of output file (should end in .geo)")
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.geomdef), "Geometry definition file does not exist: %s"%cmdline.geomdef
  assert osp.isfile(cmdline.paramdef), "Parameter definition file does not exist: %s"%cmdline.paramdef

  #Read in the two yaml files
  geomdef=useful.readyaml(cmdline.geomdef)
  paramdef=useful.readyaml(cmdline.paramdef)
  
  #Create the file
  write_one_geo(geomdef, paramdef, cmdline.geofile)
