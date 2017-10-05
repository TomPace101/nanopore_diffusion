#Process the jinja2 template into one (or more) gmsh .geo file(s)
#Usage:
#python buildgeom.py geomdef paramdef
#for more details, see
#python buildgeom.py -h

## TODO: validation of geometric inputs (different formulas for different geometries)

#Standard library
import argparse
import os.path as osp
import sys

#Site packages
from jinja2 import Environment, FileSystemLoader

#Local
sys.path.append(osp.abspath('..'))
import useful
from folderstructure import *

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

  #Dictionary of surface loops (and volumes)
  t_input['surfloops']=dict([(n, ', '.join(['%d'%x for x in surfs])) for n,surfs in geom.surfloops.items()])

  return t_input

def write_one_geo(geomdef, paramdef, geofile):
  """Generate a single geo file based on a geometry defintion dictionary and parameter dictionary
  Inputs:
    geomdef = geometry definition dictionary, which must contain:
      tmplfile: geometry template file
      ptdict: dictionary of points and their corresponding mesh density parameter name
      geomtable: mapping of surfaces to points
      internalsurfs: list of surfaces to exclude from surface loops
      revsurfs: list of surfaces needing orientation reversal
      nonplanar: list of surfaces that are not planar surfaces
    paramdef = parameter defintion dictionary, which must contain all parameters needed by the geometry template file
    geofile = path to output .geo file, as string
  No return value. The .geo file is written."""

  #Namepsace the geometry definition for convenience
  geom=argparse.Namespace(**geomdef)

  #Get the input dictionary for the template
  t_input = prepare_template_input(geom, paramdef)

  #Load template
  env=Environment(loader=FileSystemLoader(['.','./mesh']), trim_blocks=True)
  tmpl=env.get_template(geom.tmplfile)

  #Render template
  outdat = tmpl.render(t_input)

  #Output result
  with open(geofile,'w') as fp:
    fp.write(outdat)
  return

def process_mesh_params(params):
  """Turn a mesh control parameters object into the arguments needed for write_one_geo and doit file dependencies
  Inputs:
    params = an object (such as a Namespace) with at least the following attributes:
    meshname = stem name for the .geo, .msh and .xml files (this script only generates the .geo file)
    lattice = stem name of the geometry defintion yaml file (such as body-centered.yaml or face-centered.yaml),
      this file is specific to the jinja2 template used to generate the .geo file
    and all the geometric parameters needed by the jinja2 template to generate the .geo file.
  Returns:
    geomyaml = geometry definition yaml file path, as string
    geomdef = geometry definition dictionary
    geofile = output .geo file path, as string"""
  geomyaml=osp.join(meshfolder,params.lattice+'.yaml') #geometry definition yaml file
  geomdef=useful.readyaml(geomyaml) #geometry definition dictionary
  geofile=osp.join(geofolder,params.meshname+'.geo') #.geo file
  return geomyaml, geomdef, geofile

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  mesh_params_docstring="""parameter definition file for the mesh
  This is a potentially multi-doc yaml file, where each document specifies one mesh to generate.
  The required attributes are as specified in the argument of process_mesh_params"""
  parser = argparse.ArgumentParser(description='Create gmsh .geo file(s)')
  parser.add_argument('mesh_params_yaml', help=mesh_params_docstring)
  cmdline=parser.parse_args()
  assert osp.isfile(cmdline.mesh_params_yaml), "Mesh parameter definition file does not exist: %s"%cmdline.mesh_params_yaml

  #Read in the yaml file
  meshruns=useful.readyaml_multidoc(cmdline.mesh_params_yaml)
  
  #Generate each requested .geo file
  for run in meshruns:
    params=argparse.Namespace(**run)
    geomyaml, geomdef, geofile = process_mesh_params(params)
    print(geofile)
    write_one_geo(geomdef, run, geofile)
