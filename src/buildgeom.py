#Process a jinja2 template into one (or more) gmsh .geo file(s)

#Standard library
from __future__ import print_function, division #Python 2 compatibility
import os
import os.path as osp
import sys

#Site packages
from jinja2 import Environment, FileSystemLoader

#Local
import folderstructure as FS
import common

#Path to this code file (for dependency list)
thisfile=sys.modules[__name__].__file__

class GeometryDefinition(common.ParameterSet):
  """Subclass of common.ParameterSet to store the problem geometry without reference to physical dimensions
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
  _required_attrs=list(__slots__)
  _folders={'tmplfile':FS.geotemplates_folder}
  _inputfile_attrs=['sourcefile','tmplfile']
  _more_inputfiles=[osp.join(FS.geotemplates_folder,'common.geo.jinja2')]

class MeshParameters(common.ParameterSet):
  """Subclass of common.ParameterSet to store the data for generating a mesh in gmsh of the problem geometry
  These parameter files are usually stored in the location specified by folderstructure.params_mesh_folder.
  Attributes:
    To be read in:
      meshname = stem name for the .geo, .msh and .xml files
      geomdefname = stem name of the geometry defintion yaml file (such as body-centered.yaml or face-centered.yaml)
        The geometry definition yaml file stores the attributes of a GeometryDefinition instance
        This file is specific to the jinja2 template used to generate the .geo file.
      tmplvalues = dictionary of physical dimensions, mesh density parameters, etc
        Note that the geometry definition file contains a list of the variables that must be defined in this dictionary.
    To be generated by methods:
      geomdef = instance of GeometryDefinition corresponding to geomdefname
      geofile = name of .geo output file, as string (not full path)"""
  __slots__=('meshname','geomdefname','tmplvalues','geomdef','geofile','_folders')
  _required_attrs=['meshname','geomdefname','tmplvalues']
  _config_attrs=['basename']+_required_attrs
  #don't need sourcefile as input file due to config
  _more_inputfiles=[thisfile,common.__file__]
  _outputfile_attrs=['geofile']
  _child_attrs=['geomdef']
  _taskname_src_attr='meshname'
  #Remember loaded geometry definitions so we aren't reading the same files over and over when generating multiple meshes
  loaded_geomdefs={}

  def __init__(self,**kwd):
    #Initialization from base class
    super(MeshParameters, self).__init__(**kwd)
    #Get folders
    self._folders={'geomdefname':FS.geomdef_folder,'geofile':osp.join(FS.geofolder,self.basename)}
    #Get name of output file
    self.geofile=self.meshname+'.geo'
    #Load the geometry definition, unless already loaded
    if not self.geomdefname in self.loaded_geomdefs.keys():
      self.loaded_geomdefs[self.geomdefname]=GeometryDefinition.from_yaml(self.full_path('geomdefname')+'.yaml')
    self.geomdef=self.loaded_geomdefs[self.geomdefname]

  def run(self):
    print(self.geofile)
    write_one_geo(self.geomdef,self,self.full_path('geofile'))

#From mapping of surfaces to points, generate:
# - mapping of loops to line and circle names
# - mapping of line and circle names to list of points
def add_entity(tup, tdict, looplist, nameprefix):
  """Create lines/circles needed for line loops, unless they already exist.
  Inputs:
    tup = line/circle tuple of point numbers (as integers) (lines have 2 points, circles have 3)
    tdict = dictionary storing the line/circle tuples
    looplist = list of lines/circles in this line loop
    nameprefix = eg 'C' for circle, 'L' for line
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
        add_entity(ctup,circles,loops[surfnum],'C')
      else:
        ltup=(startpt,pttup[indx])
        add_entity(ltup,lines,loops[surfnum],'L')
      #Next point
      startpt=pttup[indx]
      indx += 1

  #Generate physical lines, if needed
  if geom.dimensions == 2:
    t_input['phys_lines']={}
    for lname in lines.keys():
      pts=lname[1:].split('_')
      pts.sort()
      lnum=''.join(pts)
      t_input['phys_lines'][lnum]=lname

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
  env=Environment(loader=FileSystemLoader([geomdef._folders['tmplfile'],'.']), trim_blocks=True)
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

#Support command-line arguments
if __name__ == '__main__':
  program_description='Create gmsh .geo file(s)'
  input_file_description="""Path to parameter definition file for the mesh
    This is a potentially multi-doc yaml file, where each document specifies one mesh to generate.
    Each document must provide the attributes for an instance of the MeshParameters class,
    or a subclass thereof appropriate for the specified geometry defintion."""
  
  common.run_cmd_line(program_description,input_file_description,MeshParameters)
  