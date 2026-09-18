This directory has a few programs related to resolution.  Simul2017 
produces full resolution files only when ires>0.  Use ires=2 to compute 
only on the last iteration.  Full resolution adds CPUtime and makes a 
big file.

DWS is always computed for the inversion nodes. The most useful DWS is 
computed when nitmax=1. This creates a full DWS file for all nodes even 
those that are fixed, so it is quick and easy to compute for evaluating 
data distribution.

When resolution is selected: nodes, resolution and DRE (diagonal 
element) files are created.  The nodes file points from the inversion 
parameters (that are described by the resolution) to the model 
parameters.  The inversion only includes non-fixed, non-linked nodes 
and also excludes nodes without adequate observations. The model 
generally has many more parameters.

reswid4 is a program that computes the spread function (also called 
resolving width by Michelini).  resconty and rescontz compute the 
smearing contours, as shown in Reyners et al. (1999).

drelink is a program that assigns the diagonal resolution element of 
each master node to all associated linked nodes.  This is useful for 
plotting dre when linked nodes are used to have larger areas related to 
some master nodes.

