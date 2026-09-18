simul2023 is an update of simul2017, January 2024

       The program 'simul2023' provided in this zenodo site is a flexible 
program for simultaneous inversion of travel-time data for 3-D velocity and 
hypocenters.  We have used it extensively since the 1980s, and it has evolved 
over that time with modifications to the original code of Thurber (1983).  
Eberhart-Phillips (1990) developed the Simulps code to use P and S travel-
times to solve for Vp and Vs.  From Simulps, other groups have evolved simul 
versions with useful modifications (Haslinger and Kissling, 2001), which are 
not included in simul2023.  All S travel-time ray-paths are calculated using 
the 3-D Vs structure in the Eberhart-Phillips (1990) code and all subsequent 
modified codes such as Eberhart-Phillips and Reyners (1997, 2012) and 
Eberhart-Phillips and Fry (2017).  
       
       Simul2023 is an update of simul2017.  This has improvements to shallow 
derivatives for shorter period group velocity, earthquake residual weighting 
and earth-flattening transformation velocity partial derivatives. An additional 
hypocentre output file includes many error parameters, which were being 
computed but not retained, and assigns quality codes (Eberhart-Phillips and 
Reyners, 2023). Hypocenter depths may be fixed on relocation, if desired.  
For consistency with simul2000, control parameters kspr and ieft may be adjusted.

       These parameterize velocity on a 3-D grid of nodes with velocity 
linearly interpolated between nodes, allowing flexible gridding through 
linking of nodes (Thurber and Eberhart-Phillips, 1999).  The solution is 
obtained by damped least squares, with no other smoothing applied, to obtain 
a model that has few artefacts and stays close to the initial model where 
there is low resolution.  For typical earthquake travel-time data, the Vs 
model is poorly constrained relative to the Vp model, less representative of 
crustal structure and difficult to use for interpreting Vp/Vs (Eberhart-
Phillips, 1989).  Thus, as described by Eberhart-Phillips and Reyners (1997), 
it is preferable to solve for Vp and Vp/Vs, when using local earthquake 
travel-time data.  This parameterization is retained for group velocity, with 
the Herrmann (2013) Vp kernels related to Vp model inversion parameters and 
Vs kernels related to partial derivatives of Vp and Vp/Vs model parameters 
(Eberhart-Phillips and Fry, 2017; Eberhart-Phillips et al., 2022).
       
       Rietbrock (2001) implemented Q inversion in simul, using t* 
observations which describe path attenuation from spectral decay.  This has 
been extremely useful in many settings, since Q has ability to characterize 
many features more distinctly than velocity (eg: Eberhart-Phillips et al., 
2020; Eberhart-Phillips, 2016).  It can be used for Qp or Qs, and does not 
need to solve for hypocenters. The code was extended for joint teleseimic and 
local inversion for velocity or Q (Eberhart-Phillips and Fry, 2018).
       
       The simul2023 code is a complete fortran code.  The beginning of the 
code has descriptions of input parameters.  An old users manual for an early
version may be of some use (Evans et al., 1994). Our work has taken the approach 
of adding new features, without eliminating previous options.  There are 
examples in the examples directory.  For further knowledge of input formats, 
the input subroutines can be reviewed. The include file has parameters 
related to array size and the user should consider adjusting these. There 
should be no need to edit the fortran code.  A notable input parameter is 
'nitmax', specifying the number of iterations, relocation only (0), or create 
synthetic data (-1).  The coordinate system model cartesian coordinates and 
distances are computed with Transverse Mercator conversion, and earth-
flattening transformation is used. The coordinates are defined at the 
beginning of input file 2 (stations), and older options of short-distance 
conversion or NZMG49 are allowed.  The model 0-km depth equates to a present 
sea-level datum, the z= -1 km depth grid is 1 km above sea-level, and 
velocity is linearly interpolated between gridded nodes. Travel-time ray-
tracing includes station elevations.  There are several options of travel-
time data format, selected by input parameter kttfor, and used in subroutine 
input4. 
       
Eberhart-Phillips, D. (1989), Investigations of crustal structures and active 
tectonic processes in the Coast Ranges, Central California, 209 pp, 
Stanford Univ., Stanford, CA.
Eberhart-Phillips, D. (1990), Three-dimensional P- and S-velocity structure 
in the Coalinga region, California, J. Geophys. Res, 95, 15343-15363.
Eberhart-Phillips, D. (2016), Northern California seismic attenuation: 3-D Qp 
and Qs models, Bull. Seimol. Soc. Am., 106, 2558-2573, 
doi:10.1785/0120160060.
Eberhart-Phillips, D., S. Bannister, and M. Reyners (2020), Attenuation in 
the mantle wedge beneath super-volcanoes of the Taupo Volcanic Zone, New 
Zealand, Geophys. J. Int., 220, 703-723, doi:doi.org/10.1093/gji/ggz455.
Eberhart-Phillips, D., and B. Fry (2017), A new scheme for joint surface wave 
and earthquake travel-time inversion and resulting 3-D velocity model for 
the western North Island, New Zealand, Phys. Earth and Plan. Int., 269, 
98–111, doi:10.1016/j.pepi.2017.05.014.
Eberhart-Phillips, D., and B. Fry (2018), Joint local and teleseismic 
inversion for 3-D velocity and Q in New Zealand, Phys. Earth and Plan. 
Int., 283, 48-66, doi:https://doi.org/10.1016/j.pepi.2018.08.005.
Eberhart-Phillips, D., and M. Reyners (1997), Continental subduction and 
three-dimensional crustal structure:  The northern South Island, New 
Zealand, J. Geophys. Res, 102, 11,843-811,861.
Eberhart-Phillips, D., and M. Reyners (2012), Imaging the Hikurangi plate 
interface region with improved local-earthquake tomography, Geophysical 
J.Int., 190, 1221-1242, doi: 1210.1111/j.1365-1246X.2012.05553.x.
Eberhart-Phillips, D. & Reyners, M., (2023),Catalogue of 2001-2011 New Zealand 
earthquakes relocated with 3-D seismic velocity model and comparison to 2019-2020 
auto-detected earthquakes in the sparsely instrumented southern South Island, 
N Z J Geol Geophys.,66, 646-653, doi: 10.1080/00288306.2022.2089171.
Eberhart-Phillips, D., P. Upton,  D. J. A. Barrell, B. Fry, S. Bourguignon,
and E. Warren-Smith (2022), The influence of basement terranes on tectonic 
deformation: joint earthquake travel-time and ambient noise tomography of 
the southern South Island, New Zealand, Tectonics, 41, e2021TC007006, 
doi:10.1029/2021TC007006.
Evans, J., D. Eberhart-Phillips, C. H. Thurber (1994), User's manual for 
SIMULPS12 for imaging Vp and Vp/Vs:  A derivative of the "Thurber" 
tomographic inversion SIMUL3 for local earthquakes and explosions, U.S. 
Geol. Surv. Open-File Rep. 94-431, 69 pp.
Haslinger, F., and E. Kissling (2001), Investigating effects of 3-D ray 
tracing methods in local earthquake tomography, Phys. Earth Planet. 
Inter., 123, 103-114.
Rietbrock, A. (2001), P-wave attenuation structure in the fault area of the 
1995 Kobe earthquake, J. Geophys. Res., 106, 4141-4154.
Thurber, C., and D. Eberhart-Phillips (1999), Local earthquake tomography 
with flexible gridding, Computers and Geoscience, 25, 809-818.
Thurber, C. H. (1983), Earthquake locations and three-dimensional crustal 
structure in the Coyote Lake area, central California, J. Geophys. Res., 
88, 8226-8236.


