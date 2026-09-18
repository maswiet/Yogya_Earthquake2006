This directory has 5 examples that demonstrate several uses of
simul2023.  It is not designed to be exhaustive or to be a tutorial.
They are run from scripts, with multiple input files, in run directories.
The simul code allows many output files, from various researchers over
many years.
Output files that are desired to retain are moved to another directory.

In previous years, the simul program was shared by providing the code and
a couple simple examples on an ftp site.

Research results have been obtained with series of gradational inversions 
and varied linking.  Each of the examples demonstrates the inversion code 
setup for one case and some of the data, and is not intended to give a 
complete research result.

An example for Q is in directory Q_exam.  This is from Eberhart-Phillips and 
Chadwick, Raukumara, JGR 2002.  This includes resolution (ires=3).  This 
option is not generally used unless required because it greatly increases 
runtime.

An example for velocity (Vp and Vp/Vs) is in directory spi. This is taken 
from Eberhart-Phillips and Reyners, Hikurangi southern plate interface, GJI 
2012.

A large Q example is in directory ncq. This is from Eberhart-Phillips,
Northern California, BSSA 2016.  It includes autolinking regions, as also 
shown in the vlint directory.  It uses the NCalif format (kttfor=4) for 
input.

An example of joint local and teleseismic travel-time inversion is in 
directory ste.  This is from Eberhart-Phillips and Fry, NZ teleseismic, PEPI 
2018.

An example of joint eq travel-time and group velocity inversion is in 
directory otf.  This is from Eberhart-Phillips et al., Otago, Tectonics 
2022.

