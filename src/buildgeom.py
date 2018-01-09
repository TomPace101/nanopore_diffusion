#Process a jinja2 template into one (or more) gmsh .geo file(s)
#Usage:
#python buildgeom.py geomdef paramdef
#for more details, see
#python buildgeom.py -h
#Can also be used as an imported module

## TODO: validation of geometric inputs (different formulas for different geometries)

#Standard library
import argparse
import os
import os.path as osp

#Site packages
from jinja2 import Environment, FileSystemLoader

#Local
from folderstructure import *
import useful

class MeshParameters(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store the data for generating a mesh in gmsh of the problem geometry
  These parameter files are usually stored in the location specified by folderstructure.params_mesh_folder.
  Attributes:
    meshname = stem name for the .geo, .msh and .xml files
    geomdef = stem name of the geometry defintion yaml file (such as body-centered.yaml or face-centered.yaml)
      The geometry definition yaml file stores the attributes of a GeometryDefinition instance
      This file is specific to the jinja2 template used to generate the .geo file.
    tmplvalues = dictionary of physical dimensions, mesh density parameters, etc
      Note that the geometry definition file contains a list of the variables that must be defined in this dictionary."""
  __slots__=('meshname','geomdef','tmplvalues')

class GeometryDefinition(useful.ParameterSet):
  """Subclass of useful.ParameterSet to store the problem geometry without reference to physical dimensions
  These parameter files are usually stored in the location specified by folderstructure.geomdef_folder.
  Attributes:
    dimensions = number of spatial dimensions (i.e. 2 for 2D mesh, 3 for 3D mesh)
    tmplfile = geometry template file, usually stored in the location specified by folderstructure.geotemplates_folder.
    tmplvars = mesh parameter variables needed by the geometry template file
    ptdict = dictionary of points and their corresponding mesh density parameter name
    geomtable = mapping of surfaces to sequence points
    surfloops = mapping of surface loops to sequence of surfaces
    nonplanar = list of surfaces that are not planar surfaces"""
  __slots__=('dimensions','tmplfile','tmplvars','ptdict','geomtable','surfloops','nonplanar')

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
  Inputs:
    geom = instance of GeometryDefinition
    paramdef = instance of MeshParameters
  Returns:
    t_input = the input dictionary for the template"""
    
  #Check that the necessary variables are defined
  missing=[var for var in geom.tmplvars if not var in paramdef.tmplvalues.keys()]
  extra=[var for var in paramdef.tmplvalues.keys() if not var in geom.tmplvars]
  assert len(missing)==0 and len(extra)==0, "Missing or extra template variables. missing=%s extra=%s"%(str(missing),str(extra))
  
  #Put geometric and mesh refinement parameters into template input
  t_input=dict(paramdef.tmplvalues)
  
  #Put dimensions into mesh
  t_input['dimensions']=geom.dimensions

  #Dictionary of points
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
  """Generate a single geo file based on a geometry defintion and parameters
  Inputs:
    geomdef = instance of GeometryDefinition
    paramdef = instance of MeshParameters
    geofile = path to output .geo file, as string
  No return value. The .geo file is written."""

  #Get the input dictionary for the template
  t_input = prepare_template_input(geomdef, paramdef)

  #Load template
  env=Environment(loader=FileSystemLoader([geotemplates_folder,'.']), trim_blocks=True)
  tmpl=env.get_template(geomdef.tmplfile)

  #Render template
  outdat = tmpl.render(t_input)

  #Create output directory if needed
  if not osp.isdir(osp.dirname(geofile)):
    os.makedirs(osp.dirname(geofile))
  #Output result
  with open(geofile,'w') as fp:
    fp.write(outdat)
  return

def process_mesh_params(params):
  """Turn a mesh control parameters object into the arguments needed for write_one_geo and doit file dependencies
  Inputs:
    params = a MeshParameters instance
  Returns:
    geomyaml = geometry definition yaml file path, as string
    geomdef = instance of GeometryDefinition
    geofile = output .geo file path, as string"""
  geomyaml=osp.join(geomdef_folder,params.geomdef+'.yaml') #geometry definition yaml file
  geomdef=GeometryDefinition.from_yaml(geomyaml)
  geofile=osp.join(geofolder,params.basename,params.meshname+'.geo') #.geo file
  return geomyaml, geomdef, geofile

def geo_from_MeshParams(params):
  """Generate the gmsh .geo file from the specified MeshParams
  Inputs:
    params = as MeshParameters instance
  No return value."""
  geomyaml, geomdef, geofile = process_mesh_params(params)
  print(geofile)
  write_one_geo(geomdef, params, geofile)

#Support command-line arguments
if __name__ == '__main__':
  #Process command-line arguments
  program_description='Create gmsh .geo file(s)'
  input_file_description="""Parameter definition file for the mesh
    This is a potentially multi-doc yaml file, where each document specifies one mesh to generate.
    Each document must provide the attributes for an instance of the MeshParameters class,
    or a subclass thereof appropriate for the specified geometry defintion."""
  
  useful.run_cmd_line(program_description,input_file_description,MeshParameters,geo_from_MeshParams)
  