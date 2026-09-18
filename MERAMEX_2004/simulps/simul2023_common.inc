c23456789a123456789b123456789c123456789d123456789e123456789f123456789g12
c   common blocks for program SIMUL2023
c     This has more variables for Rayleigh wave velocity observations
c     from surface waves 'gsw'.
c
c  definitions:
c   nsts = number of stations
c   neqs = number of earthquakes
c   nbls = number of blasts
c   nsht = number of shots
c   ntel = number of teleseismic events
c   ngsw = number of group (and/or phase) velocity Rayleigh surface wave obs
c   nobs = number of observations for a single event.
c   nx  = number of nodes in the x direction.
c   ny = number of velocity nodes in y direction.
c   nz = number of velocity nodes in z direction.
c   nodes2=(nx-2)(ny-2)(nz-2) - number of nodes which are inverted for.
c   nparv=number of velocity parameters, 2*nodes2 if both P and S
c         otherwise=nodes2
c   nparvi=number of velocity parameters allowd to invert for
c   nstsi = number of station to invert for station terms
c   npar = nparv ( station delays not included)
c         nparv + nsts  (station delays,P only)
c         nparv + 2*nsts (station delays,P and S)
c   npari = nparvi (station delays not included)
c           nparvi + nstsi (station delays, P only)
c           nparvi + 2*nstsi (station delays,P and S)
cfhdmep
c   nparpi = number of P-velocity nodes allowed to invert for (reduced by fixed nodes)
c   nparsi = number of S-velocity nodes allowed to invert for (reduced by fixed nodes)
c   nrowp,nrows,nrowst = mbla(2,3,4) for output to file 18 (full resol. matrix)
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c  TO CHANGE ARRAY SIZES,
C  PARAMETERS IN THIS SECTION CAN BE CHANGED
c
c   maxev: maximum number of events (earthquakes+blasts+shots+tele)
      parameter(maxev=2500)
c   maxtev: maximum number of teleseismic events
      parameter(maxtev=150)
c   maxgsw: maximum number of group (and/or phase) velocity Rayleigh observations
      parameter(maxgsw=1000)
      parameter(mxevtot=maxev+maxgsw)
c
c   maxobs: maximum number of arrivals for a particular event.  Note
c          that maxobs can generally be set much lower than maxsta
c          since events will not generally be recorded at all
c          stations.
      parameter(maxobs=0605)
c
c   maxclu: maximum number of clusters
      parameter(maxclu=80)
c   maxcev: maximum number of events per cluster
      parameter(maxcev=36)
      parameter(mxhparc=maxcev*4)
c   maxedt: maximum number of eq-differential times per cluster
      parameter(maxedt=1600)
c   mxobsa: maximum number of obs for eq or edt
      parameter(mxobsa=maxobs+maxedt)
c mxobsc: maximum obs for cluster, includes both ttime and 
c    differential time. Usually will not have maxobs for all
c    events in cluster.  mxobsc <= maxcev*maxobs + maxedt
      parameter(mxobsc=2000)
c
c   maxsta: maximum number of stations in station list
      parameter(maxsta=24900)
c
c   mxpari: maximum number of parameters to invert for.  This 
c           includes Vp and Vp/Vs nodes and P and S station 
c           corrections.  mxpari will generally be much lower
c           than maxpar because velocity nodes and station
c           corrections can be fixed. (Station corrections
c           are generally fixed to zero in the region of the
c           velocity inversion.)
c      parameter(mxpari=9100)
      parameter(mxpari=10950)
c      parameter(mxpari=11550)
c
c   maxpar: maximum number of potential parameters that could
c           be included in the inversion.  This includes all 
c           velocity nodes (includes both Vp and Vp/Vs) that
c           define the 3-D velocity model (nparv). Does not 
c           include edge nodes.
c           If station corrections are included in the inversion,
c           maxpar also includes stations (nparv+nsts or 
c           nparv + 2*nsts if P and S-P are used).
c           maxpar is set much larger than mxpari to allow
c           for assigning regional velocity (from prior
c           coarse inversions) values (fixed) in peripheral
c           areas where there is nil resolution.
      parameter(maxpar=250000)
c
c   maxnx: maximum number of nodes in x direction
      parameter(maxnx=80)
c   maxny: maximum number of nodes in y direction
      parameter(maxny=140)     
c   maxnz: maximum number of nodes in z direction
      parameter(maxnz=80)     
c  maxpag: maximum number of x or y values for box to compute gsw
c          velocity and partials
      parameter(maxpag=100)
c  mxnpsksw: maximum number of points for each gsw sensitivity kernel
      parameter(mxnpsk=300)
c
c
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c
c  common block variables:
c
c   control parameters
c
      common/contrl/ nsts,neqs,nobt,wnobt,nobtt,
     * ires,nitloc,eigtol,rmscut,cvpmax,cvpmin,
     * dvpmx,dvsmx,hitct,nsht,nbls,nevt,i3d,nitmax,kspr,snrmct,ihomo,
     * rmstop,delt1,delt2,ddlt,res1,res2,res3,dres12,dres23,
     * stepl,ieft,wtsht,wtsp,kout,kout2,kout3,kttfor,ifixz
     * ,ndip,iskip,scale1,scale2,iusep,iuses,invdel,zmin,dxmax,idmp
     * ,vdamp(4),ifixl,xfac,tlim,nitpb(2),nobtp,nobts,nobtex,nobteq
     * ,nwrt,nswrt,rderr,ercof,wnobtp,wnobts,iuseq,iuse2t,iqwt,qrmax,
     *  nebs,ntel,wttel,ntobmin,reste1,reste2,reste3,dreste12,dreste23,
     *  nobtte,qmin,
     *  ngsw,dhsw,boxmin,zswmax,wtgsw,wnobtg
c
c   control parameters for receiver-pair differential travel-time
      common/conrdt/ nrdtsh,iuserdt,rdisrdt,facrayr,frayrdt,gminrdt,
     *  mxiqrdt,wtrdtsht,res4,nordt,nordtp,nordts
c
c   control parameters for cluster earthquake-pair differential
c     travel-time
      common/conclu/ nclu,iuseclu,edisedt,facraye,frayedt,gminedt,
     *  mxiqedt,wtepdt,noedt,noedtp,noedts,nobtce,nstaepair,wrmscemx,
     * gminedt0,minedt
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c paramaters to define hypocenter quality codes, subroutine qualhy
c put in include so can easily be changed
      parameter(serzfix=20.0)
      character*1 qcode(5)
      data qcode/'A','B','C','D','E'/
      integer nobqc(5),nsqc(5),igapqc(5),idltqc(5)
      real rmsqc(5),rzdmqc(5),erhqc(5),erzqc(5)
      data nobqc/15,12,8,6,4/
      data nsqc/4,2,1,0,0/
      data igapqc/140,180,240,270,360/
      data rmsqc/0.2,0.3,0.4,0.6,10./
      data rzdmqc/0.5,0.3,0.1,0.0,0.0/
      data erhqc/1.5,2.5,4.,8.,200./
      data erzqc/2.5,4.,6.,12.,200./
c set all idltqc big, since not using, and relying on rzdm
      data idltqc/500.,500.,500.,500.,500./
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   event observation and station arrays
      character*1 ecns,ecew,scns,scew,tecns,tecew,cqhy
      character*11 evid
      common/events/ mino(maxev),seco(maxev),ltde(maxev),eltm(maxev),
     *  ihr(maxev),lnde(maxev),elnm(maxev),evc(3,maxev),kobs(maxev),
     *  iyrmo(maxev),iday(maxev),rmag(maxev),wrmsr(maxev),seco2(maxev),
     * igap(maxev),idltmn(maxev),mina(maxev),seca(maxev),kobps(2,maxev),
     *  elat(maxev),elon(maxev),ecns(maxev),ecew(maxev),evid(maxev),
     *  serhy(5,maxev),iyr2(maxev),imo(maxev),cqhy(maxev),
     *  iiday(maxev),iihr(maxev)
      common/ntemp/ netemp,nbtemp,neb
      common/sepdat/ totrms(mxevtot),seprms(mxevtot)
      common/fastrc/ jfl,nco(maxobs,maxev),ndo(maxobs,maxev)
c
      character*5 stn5
       character*2 net
      character*6 stn6
      common/observ/ stn(maxsta),ltds(maxsta),sltm(maxsta),lnds(maxsta),
     *  slnm(maxsta),stc(3,maxsta),isto(maxobs,maxev),wt(maxobs,maxev),
     *  secp(maxobs,maxev),intsp(maxobs,maxev),pdl(maxsta),sdl(maxsta),
     *  dlta(maxobs,maxev),nfixst(maxsta),iclock(maxsta),
     *  rdlta(maxobs,maxev),stn6(maxsta),slat(maxsta),slon(maxsta),
     *  scns(maxsta),scew(maxsta),stn5(maxsta),net(maxsta)
      common/obsiw/ iw(maxobs,maxev)
      common/nrdsta/nrd(maxsta,2)
      character*4 rmk,rmkc,rmkrev
      common/obsrmk/ rmk(maxobs,maxev),rmkc(maxobs,maxcev,maxclu),
     2  rmkrev(maxobs,maxev)
      common/fpsdat/ az(maxobs),toa(maxobs)
c
c
c teleseismic event and observ arrays
      common/obstel/tepp(3,maxobs,maxtev),tecen(3,maxtev),
     *  teplat(maxobs,maxtev),teplon(maxobs,maxtev),
     *  tepltm(maxobs,maxtev),teplnm(maxobs,maxtev),telpad(maxtev),
     *  telpads(maxtev),secte(maxobs,maxtev),tecnlt(maxtev),
     *  tecnln(maxtev),ltdtep(maxobs,maxtev),lndtep(maxobs,maxtev),
     *  tecns(maxobs,maxtev),tecew(maxobs,maxtev)
c
c   receiver-pair differential travel-time arrays
      common/obsrdt/ rdtsec(maxobs,maxev), intsprd(maxobs,maxev),
     *  istord(2,maxobs,maxev), wtrdtobs(maxobs,maxev),
     *  iwrdt(maxobs,maxev), jobsrd(2,maxobs,maxev),kobsrdt(maxev),
     *  nrprdt(2,maxobs,maxev), wrmsrdt(maxev),resrdt(maxobs),
     *  ttmrdt(maxobs),wtcombrd(maxobs,maxev)
c
c   cluster earthquake-pair dt event and observ arrays
      common/cevents/ minoce(maxcev,maxclu),secoce(maxcev,maxclu),
     *  ltdce(maxcev,maxclu),celtm(maxcev,maxclu),ihrce(maxcev,maxclu),
     *  lndce(maxcev,maxclu),celnm(maxcev,maxclu),cevc(3,maxcev,maxclu),
     *  kobsce(maxcev,maxclu),iyrmoce(maxcev,maxclu),ncev(maxclu),
     *  ccenc(3,maxclu),idayce(maxcev,maxclu),rmagce(maxcev,maxclu),
     *  wrmsrce(maxcev,maxclu),igapce(maxcev,maxclu),
     *  idltmnce(maxcev,maxclu),minace(maxcev,maxclu),
     *  secace(maxcev,maxclu),nobsc(maxclu),nobsct(maxclu),
     *  serhyc(4,maxcev,maxclu)
      common/cobserv/ istoc(maxobs,maxcev,maxclu),
     *  wtc(maxobs,maxcev,maxclu),secpc(maxobs,maxcev,maxclu),
     *  intspc(maxobs,maxcev,maxclu),iwc(maxobs,maxcev,maxclu),
     *  dltac(maxobs,maxcev,maxclu),rdltac(maxobs,maxcev,maxclu)
      common/cobsedt/edtsec(maxedt,maxclu),intsped(maxedt,maxclu),
     *  wtedtobs(maxedt,maxclu),istoed(2,maxedt,maxclu),
     *  jeved(2,maxedt,maxclu),
     *  iwedt(maxedt,maxclu),jobsed(2,maxedt,maxclu),
     *  kobsedt(maxclu),nrpedt(2,maxedt,maxclu),
     *  wrmsedt(maxclu),wrmsclu(maxclu),
     *  resedt(maxedt),ttmedt(maxedt),wtcombed(maxedt,maxclu)
      common/cfastrc/ jflc(maxcev),ncoc(maxobs,maxcev,maxclu),
     *  ndoc(maxobs,maxcev,maxclu)
c
c  group (and phase) velocity surface wave sensitivity kernel arrays
      common/skerngsw/ sksw(mxnpsk,2),zsksw(mxnpsk),
     2 nsk,persksw,nsku,wsksw(mxnpsk,2)
c  group velocity surface wave observation arrays
      common/obsgsw/ vgsw(maxgsw),pergsw(maxgsw),iqgsw(maxgsw),
     2 wtswobs(maxgsw),rlatsw(maxgsw),rlonsw(maxgsw),elevsw(maxgsw),
     3 evsw(3,maxgsw),xpag(maxpag,maxgsw),ypag(maxpag,maxgsw),
     4 nxpag(maxgsw),nypag(maxgsw),wxpag(maxpag,maxgsw),icgsw(maxgsw),
     5 wypag(maxpag,maxgsw),vgswwm(maxgsw),
     6 wrmsgsw(maxgsw),vgswm(maxgsw),rsvgsw(maxgsw)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   velocity arrays
c
c   /vmod3d/
c    vel(nx,ny,2*nz)     ! 1st half is P velocity, 2nd half is S velocity
      parameter(maxnz2=maxnz*2)
c
      common/vmod3d/bld,nx,ny,nz,xn(maxnx),yn(maxny),zn(maxnz),
     *  vel(maxnx,maxny,maxnz2),vpvs(maxnx,maxny,maxnz),
     *  qval(maxnx,maxny,maxnz)
c
c   /locate/
c   ixkms: size of map in x direction.  The distance between 
c            furthest nodes.  This will be in units of km if bld=1.0 
c            or in units of 0.1km if bld=0.1
      parameter(ixkms=10000)
c   iykms=size of map in y dir.
      parameter(iykms=10000)
c   izkms=size of map in z dir.
      parameter(izkms=10000)
c
      common/locate/ xl,yl,zl,ixloc(ixkms),iyloc(iykms),izloc(izkms)
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  inversion arrays
c
      common/hypinv/ res(mxobsa),dth(maxobs,5),dthp(maxobs,5),
     *  wtcomb(maxobs,maxev),ttm(maxobs)
      common/chypinv/ resc(maxobs,maxcev),dthc(maxobs,4,maxcev),
     *  wtcombc(maxobs,maxcev,maxclu),ttmc(maxobs,maxcev)
c
cfhdmep nparpi, nparsi, nrowp,nrows,nrowst added here
c   /modinv/
c      parameter(mxdtm=maxobs*mxpari)
      parameter(mxdtm=mxobsa*mxpari)
      common/modinv/ dtm(mxdtm),dtmp(mxdtm),resp(mxobsa),nz2,nx1,
     *  ny1,nz1,nparvi,nstsi,npari,inf,nfix(maxpar),ndexfx(maxpar),
     *  mdexfx(maxpar),nparv,npar,nodes,nxy,hit(maxpar),hitall(maxpar),
     *  vadj(maxpar),khit(maxpar),nx2,nxy2,nodes2,sumhit(maxnz2),
     *  qadj(maxpar),nparpi,nparsi,nrowp,nrows,nrowst,nparvs,nrowte,
     *  imerge(maxpar),jequal(maxpar),ltype(maxpar),infl,imast,ilink,
     *  nparhy,nparvsi
      character*1 cnode,canode
      common/modinc/cnode(maxpar),canode(maxpar)
c
c   /solutn/
      parameter(mgsol=(mxpari*(mxpari+1))/2)
      parameter(mxpri1=mxpari+1)
      common/solutn/ g(mgsol),g1(mgsol),rhs(mxpri1),rhs1(mxpri1),
     *  index(maxpar),jndex(maxpar),mbl,nbtot,ssqr,mbl1,ssqr1,
     *  ssqrw,ssqrw1,ssqrwp,ssqrws,ssqrwg,var,varw,var1,varw1,
     *  ndof,ndof1,wndof,wndof1,sumnrm(5)
c
      common/resltn/ covdi(maxpar),drm(maxpar),stderr(maxpar)
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c
c  DIMENSIONS FOR RAY TRACING
c   
c    NOTE: This is not designed to be changed by the user.
c          There are other dimension and commons throughout
c          the program besides /rpath/.  /rpath/ is just
c          included here because it also uses maxobs to
c          dimension the array.
c
c    nseg = number of segments in longest ray. was 256 segments (23jan06)
c
c     nseg=2**(1+nint(3.32193*alog(sep/scale1))). The
c     maximum number of segments allowed (see subroutine setup)
c     is now 2**n2exp, where n2exp is below in include.
c
c     ncr=1+.5*sep/scale2   - - the maximum number of curves
c                        for the circular rays.
c
c     sep=max dist between event and station. Scale1 is set in
c     input parameters. So is scale2.
c
c    (3*nseg) : strpth,fstpth,trpth1
c    (nseg) : pthsep
c    (3,9) : dipvec
c    (3*nseg,9) : disvec,trpath
c    (9) : trtime
c   /rpath/
c   (3,nseg,nobs) : rp
c   (nobs) : nrp,pl,tt
c
      parameter(n2exp=8)
      common/rpath/rp(3,260,maxobs),nrp(maxobs),ttc(maxobs),
     *  xmov,pl(maxobs)
c
c also save ray paths for obs of events in cluster
      common/crpath/ rpce(3,260,maxobs,maxcev),nrpce(maxobs,maxcev),
     *  ttcce(maxobs,maxcev),plce(maxobs,maxcev)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
