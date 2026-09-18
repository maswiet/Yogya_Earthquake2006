This directory has codes for setting up the input velocity model.
NOTE that we edit the velocity input files, to add in various fixed nodes
and linked nodes.  Also the beginning of the velocity input file
can have comment lines, with 'c' as the first character.
These various programs are useful to manipulate velocity files and sections.
They are not necessarily written for usefulness to other people,
so there are no directions or guarantees.

'vlinterp' is a program to interpolate velocity grids.
'vlinterp2' is another that uses different orgin and rotation.
It is an old program that does not have general TM coordinates
and only has NZ or short-distance conversion option.
See interp2 directory.

vlavg averages together velocity files, generally used for varied
linking runs.

vlmerge merges a selected area or depth from one velocity file
into another.

Simul does an actual inversion rather than the less accurate LSQR
of many codes, and thus it is best to try to fix nodes that are 
not needed.

All programs are fortran, have been compiled with f77
  eg: f77 -C -O2 autolink9dxdy.f -o autolink9dxdy
Then the executable put in a bin directory which is in the person's path.

Before deciding on the areas to be linked or fixed, we generally run
some 1 iteration inversions with all nodes free to produce a 
complete DWS file.  The linking routines that we have used are based
on that, although other methods would be possible.

The approach is to have a particular inversion study or inversion 
step within a larger regional model.  That allows raypaths outside 
of the solution area, to have more stable hypoenters and more 
reliable velocity.

The area (or depths) outside the inversion area can be assigned 
fixed nodes.  Non-fixed nodes are called free nodes.
The programs vlfixpart and vlfreepart create these input sections.

Many researchers in tomography have strategies to accommodate 
uneven data distribution with uneven grids or varied size areas
related to solution parameters.  Both approaches can be done in 
simul2017 with linked nodes.

Linked nodes can have constant velocity perturbation, or
can have linear.  For the linear, their velocity has linear 
gradients between their master node and surrounding master and
free nodes.

'autolink9dxdy' links areas of size up to 3dx by 3dy, for input dx and dy.
It has been used in many studies.  Here there is an example taken
from Northern California Q (BSSA2016), in ncq_link directory.
This also used 'dwsunlink' to unlink high DWS nodes.

'link4' and 'link4diag' are programs that like 4 nodes together 
either in xy box or diagonal box.  Results are produced for 
the 4 possible patterns. They do not consider DWS.
But after running link4, the results can be adjusted with dwsunlink.

'autolink12y' is an older program that links up to 12 nodes
together, 2x, 3y, 1-2 z.  It just does one realization based on 
the starting node indices given.

'vlfixlnk' is program that can be used to correct velocity files
that have nodes assigned fixed multiple times or fixed, and also
included in linked.  This is handy if run linking and then decide
to fix parts, or run vlfreepart a few different ways for fixing. 
We like to carefully consider free, linked, and fixed nodes 
during stages of gradational inversions, maintaining the regional
3D model. In this process we like to manually construct the 
velocity input file.
For instance a linking code like link4diag might be run 
to create linked nodes, then vlfixpart might be run to create 
fixed node areas for a couple block regions that might overlap. 
Then these separate fixed and linked file sections could be 
combined together into a velocity input file that might have
nodes assigned multiple times which needs correction to run. 
For that, vlfixlnk would correct it.
