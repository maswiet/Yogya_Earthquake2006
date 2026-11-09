#!/usr/bin/env python3
"""Statistical stress inversion from focal mechanisms (Michael 1984, 1987) with
bootstrap confidence and random nodal-plane selection (as in SATSI/MSATSI).

Reads eqt/mt/results/mt_catalog.csv, selects quality mechanisms, inverts for the
reduced (deviatoric) stress tensor, derives principal stress axes (sigma1>=sigma2>=sigma3,
sigma1 = maximum compression), stress ratio R=(s1-s2)/(s1-s3), and SHmax azimuth
(Lund & Townend 2007). Bootstrap gives orientation confidence clouds.

Writes:
  eqt/mt/results/stress_result.json   (axes, R, SHmax + bootstrap percentiles)
  eqt/figures/stress_stereonet.png    (statistical stress-direction stereonet)

Pure numpy. Convention: x=North, y=East, z=Down (Aki & Richards).

Self-test:  python stress_inversion.py --selftest
Run:        python stress_inversion.py [--vrmin 0.4] [--dcmin 50] [--nmin 6] [--nboot 2000]
"""
import sys, os, json
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CAT  = os.path.join(ROOT, "mt/results/mt_catalog.csv")
OUTJ = os.path.join(ROOT, "mt/results/stress_result.json")
OUTF = os.path.join(ROOT, "figures/stress_stereonet.png")

D2R = np.pi/180.0

# ---------------------------------------------------------------- geometry
def sdr_to_normal_slip(strike, dip, rake):
    """strike/dip/rake (deg) -> (normal, slip) unit vectors, NED (x=N,y=E,z=Down).
    Aki & Richards (1980) convention. normal points to the hanging wall (up-ish)."""
    s, d, r = strike*D2R, dip*D2R, rake*D2R
    # fault normal (upward normal of the footwall / plane)
    n = np.array([-np.sin(d)*np.sin(s),
                   np.sin(d)*np.cos(s),
                  -np.cos(d)])
    # slip of hanging wall relative to footwall
    u = np.array([ np.cos(r)*np.cos(s) + np.sin(r)*np.cos(d)*np.sin(s),
                   np.cos(r)*np.sin(s) - np.sin(r)*np.cos(d)*np.cos(s),
                  -np.sin(r)*np.sin(d)])
    return n/np.linalg.norm(n), u/np.linalg.norm(u)

def aux_plane(strike, dip, rake):
    """Return the auxiliary nodal plane (strike,dip,rake) for a given plane."""
    n, u = sdr_to_normal_slip(strike, dip, rake)
    # auxiliary plane: normal = slip, slip = normal
    return nu_to_sdr(u, n)

def nu_to_sdr(n, u):
    """Inverse: normal & slip unit vectors (NED) -> strike,dip,rake (deg)."""
    # ensure normal points upward (z down -> negative z is up); flip both if needed
    if n[2] > 0:
        n = -n; u = -u
    dip = np.arccos(-n[2])
    # strike direction
    strike = np.arctan2(-n[0], n[1])
    sind = np.sin(dip)
    if abs(sind) < 1e-8:  # horizontal plane, strike undefined -> 0
        strike = 0.0
        sv = np.array([np.cos(strike), np.sin(strike), 0.0])
    else:
        sv = np.array([-n[1], n[0], 0.0]); sv /= np.linalg.norm(sv)
    # rake = angle of slip from strike vector, measured in plane
    dv = np.cross(n, sv)              # down-dip vector
    rake = np.arctan2(np.dot(u, dv), np.dot(u, sv))
    st = (np.degrees(strike)) % 360.0
    return st, np.degrees(dip), np.degrees(rake)

def axis_trend_plunge(v):
    """unit vector NED -> (trend deg from North CW, plunge deg below horizontal)."""
    v = v/np.linalg.norm(v)
    if v[2] < 0:   # force downward (lower hemisphere)
        v = -v
    plunge = np.degrees(np.arcsin(v[2]))
    trend  = np.degrees(np.arctan2(v[1], v[0])) % 360.0
    return trend, plunge

def tp_to_vec(trend, plunge):
    """(trend,plunge deg) -> unit vector NED (downward hemisphere)."""
    t, p = np.radians(trend), np.radians(plunge)
    return np.array([np.cos(p)*np.cos(t), np.cos(p)*np.sin(t), np.sin(p)])

def axis_cone(tps, mean_vec, q=68):
    """Angular radius (deg) of the smallest cone about mean_vec holding q% of
    bootstrap axes. Treats axes as undirected (antipodal-safe)."""
    m = mean_vec/np.linalg.norm(m) if (m:=mean_vec) is not None else None
    ang=[]
    for tr,pl in tps:
        v = tp_to_vec(tr,pl)
        c = abs(np.clip(np.dot(v, m), -1, 1))      # undirected
        ang.append(np.degrees(np.arccos(c)))
    return float(np.percentile(ang, q))

def mean_axis(tps):
    """Principal direction of a set of undirected axes (eigenvector of sum v v^T)."""
    M=np.zeros((3,3))
    for tr,pl in tps:
        v=tp_to_vec(tr,pl); M+=np.outer(v,v)
    w,V=np.linalg.eigh(M)
    return V[:,-1]

def circ_axial(az_deg, q_lo=16, q_hi=84):
    """Mean & 68% spread for axial data (period 180 deg), returned as
    (mean, half_width_deg). Uses doubled-angle vector averaging."""
    a=np.radians(np.asarray(az_deg)*2.0)
    C,Sn=np.cos(a).mean(), np.sin(a).mean()
    mean=(np.degrees(np.arctan2(Sn,C))/2.0)%180.0
    Rlen=np.hypot(C,Sn)
    # circular standard deviation (deg, axial)
    sd=np.degrees(np.sqrt(-2*np.log(max(Rlen,1e-9))))/2.0
    return float(mean), float(sd)

# ------------------------------------------------- Michael linear inversion
# deviatoric symmetric basis tensors (trace 0), NED
_B = [
    np.array([[1,0,0],[0,0,0],[0,0,-1]], float),
    np.array([[0,0,0],[0,1,0],[0,0,-1]], float),
    np.array([[0,1,0],[1,0,0],[0,0,0]], float),
    np.array([[0,0,1],[0,0,0],[1,0,0]], float),
    np.array([[0,0,0],[0,0,1],[0,1,0]], float),
]

def _design_rows(n):
    """3x5 design block: predicted shear-traction components for each basis tensor."""
    A = np.zeros((3,5))
    for k,B in enumerate(_B):
        t = B @ n
        tn = n @ t
        shear = t - tn*n           # tangential traction (Wallace-Bott)
        A[:,k] = shear
    return A

def invert_stress(normals, slips):
    """Least-squares reduced-stress tensor from fault normals & unit slip vectors."""
    A = np.vstack([_design_rows(n) for n in normals])
    d = np.hstack([s for s in slips])
    m, *_ = np.linalg.lstsq(A, d, rcond=None)
    T = sum(mk*B for mk,B in zip(m, _B))
    # misfit (deg): angle between predicted shear and observed slip
    ang = []
    for n,s in zip(normals, slips):
        t = T@n; sh = t - (n@t)*n
        if np.linalg.norm(sh) > 0:
            c = np.clip(np.dot(sh/np.linalg.norm(sh), s), -1, 1)
            ang.append(np.degrees(np.arccos(c)))
    return T, (np.mean(ang) if ang else np.nan)

def principal_axes(T):
    """Eigen-decompose. Returns sigma1,2,3 unit vectors (NED) and R.
    Convention: compression positive. In the traction formulation T n gives the
    traction; the most COMPRESSIVE principal stress = eigenvector of the most
    NEGATIVE eigenvalue (tractions point toward the plane under compression).
    Validated by --selftest."""
    w, V = np.linalg.eigh(T)          # ascending eigenvalues
    # w[0]<=w[1]<=w[2]; sigma1 (max compression) <-> most negative eigenvalue
    s1 = V[:,0]; s2 = V[:,1]; s3 = V[:,2]
    l1, l2, l3 = w[0], w[1], w[2]     # l1 most compressive
    R = (l1 - l2)/(l1 - l3) if (l1-l3)!=0 else np.nan
    return s1, s2, s3, R

def shmax_azimuth(T):
    """SHmax azimuth (deg from North) via Lund & Townend (2007), NED tensor."""
    # horizontal 2x2 block in (N,E); SHmax is eigenvector of max compression there,
    # but Lund & Townend give full-tensor-consistent horizontal projection.
    s1v, s2v, s3v, R = principal_axes(T)
    # Lund & Townend: build from principal axes & R
    # Represent sigma1..3 magnitudes as 1, 1-R, 0 (compression positive)
    mags = np.array([1.0, 1.0-R, 0.0])
    axes = np.array([s1v, s2v, s3v])            # rows = axes (NED)
    # reconstruct normalized stress tensor (compression positive)
    S = sum(mags[i]*np.outer(axes[i], axes[i]) for i in range(3))
    # horizontal block in N,E (indices 0,1)
    Sh = S[np.ix_([0,1],[0,1])]
    wh, Vh = np.linalg.eigh(Sh)                 # ascending
    vmax = Vh[:, -1]                            # max compression in horizontal
    az = np.degrees(np.arctan2(vmax[1], vmax[0])) % 180.0
    return az, R

# ----------------------------------------------------------------- bootstrap
def bootstrap(sdr_list, nboot=2000, rng=None):
    """Each resample: draw events with replacement, pick a random nodal plane each."""
    rng = rng or np.random.default_rng(42)
    N = len(sdr_list)
    S1=[]; S2=[]; S3=[]; RS=[]; AZ=[]
    for _ in range(nboot):
        idx = rng.integers(0, N, N)
        nn=[]; uu=[]
        for i in idx:
            st,dp,rk = sdr_list[i]
            if rng.random() < 0.5:
                st,dp,rk = aux_plane(st,dp,rk)
            n,u = sdr_to_normal_slip(st,dp,rk)
            nn.append(n); uu.append(u)
        T,_ = invert_stress(nn, uu)
        s1,s2,s3,R = principal_axes(T)
        az,_ = shmax_azimuth(T)
        S1.append(axis_trend_plunge(s1)); S2.append(axis_trend_plunge(s2))
        S3.append(axis_trend_plunge(s3)); RS.append(R); AZ.append(az)
    return (np.array(S1), np.array(S2), np.array(S3), np.array(RS), np.array(AZ))

# ----------------------------------------------------------------- self-test
def selftest():
    rng = np.random.default_rng(0)
    # impose a known stress: sigma1 (compression) N-S horizontal, sigma3 vertical
    # -> normal/thrust regime test with sigma1 horizontal N (az 0), R=0.5
    az_true = 20.0
    a = az_true*D2R
    s1_t = np.array([np.cos(a), np.sin(a), 0.0])      # horizontal, az=20
    s3_t = np.array([0,0,1.0])                        # vertical
    s2_t = np.cross(s3_t, s1_t)
    R_t = 0.5
    mags = np.array([1.0, 1.0-R_t, 0.0])              # compression positive
    axes = [s1_t, s2_t/np.linalg.norm(s2_t), s3_t]
    T_t = sum(mags[i]*np.outer(axes[i],axes[i]) for i in range(3))
    T_t = T_t - np.trace(T_t)/3*np.eye(3)             # deviatoric
    # generate faults, slip = shear-traction direction of T_t
    sdr=[]
    for _ in range(40):
        st=rng.uniform(0,360); dp=rng.uniform(25,85); rk=rng.uniform(-180,180)
        n,_u = sdr_to_normal_slip(st,dp,rk)
        t=T_t@n; sh=t-(n@t)*n
        if np.linalg.norm(sh)<1e-6: continue
        s=sh/np.linalg.norm(sh)
        st2,dp2,rk2 = nu_to_sdr(n, s)
        sdr.append((st2,dp2,rk2))
    nn=[]; uu=[]
    for st,dp,rk in sdr:
        n,u=sdr_to_normal_slip(st,dp,rk); nn.append(n); uu.append(u)
    T,mis = invert_stress(nn,uu)
    s1,s2,s3,R = principal_axes(T)
    az,_ = shmax_azimuth(T)
    t1,p1 = axis_trend_plunge(s1); t3,p3 = axis_trend_plunge(s3)
    print(f"[selftest] imposed: SHmax az={az_true}  R={R_t}")
    print(f"[selftest] recovered: sigma1 trend={t1:.1f} plunge={p1:.1f}  "
          f"sigma3 trend={t3:.1f} plunge={p3:.1f}")
    print(f"[selftest] SHmax az={az:.1f}  R={R:.2f}  mean misfit={mis:.2f} deg")
    ok_az = min(abs(az-az_true), 180-abs(az-az_true)) < 5
    ok_s1 = p1 < 10          # sigma1 should be ~horizontal
    ok_R  = abs(R-R_t) < 0.1
    print(f"[selftest] PASS az={ok_az} s1_horizontal={ok_s1} R={ok_R}")
    return ok_az and ok_s1 and ok_R

# ----------------------------------------------------------------- main
def load_quality(vrmin, dcmin, nmin, nusedmin=5):
    """Quality gate. IMPORTANT: gate on `nsta_used` (stations actually inverted),
    NOT `nsta` (stations extracted). VR anti-correlates with nsta_used (few-station
    fits are trivially high-VR but azimuthally unconstrained), so a VR-only gate
    selects unreliable mechanisms. Requires the augmented catalog (nsta_used col)."""
    import pandas as pd
    d = pd.read_csv(CAT)
    d = d[d.status=="ok"].copy()
    if "nsta_used" not in d.columns:
        raise SystemExit("catalog lacks nsta_used column; run the log-parse augmentation first")
    q = d[(d.VR>=vrmin) & (d.dc_perc>=dcmin) & (d.nsta_used>=nusedmin)].copy()
    return d, q

def plot_stereonet(q, s1,s2,s3, boot, shmax_az, R, fname):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    S1b,S2b,S3b,Rb,AZb = boot
    def tp2xy(trend, plunge):   # equal-area lower hemisphere
        tr = np.radians(trend); pl = np.radians(plunge)
        r = np.sqrt(2)*np.sin((np.pi/2-pl)/2)
        return r*np.sin(tr), r*np.cos(tr)
    fig, ax = plt.subplots(1,2, figsize=(13,6.2),
                           gridspec_kw={'width_ratios':[1.15,1]})
    a=ax[0]
    th=np.linspace(0,2*np.pi,200); a.plot(np.cos(th),np.sin(th),'k',lw=1.2)
    for t in range(0,360,30):
        x,y=np.sin(np.radians(t)),np.cos(np.radians(t)); a.plot([0,x],[0,y],color='0.9',lw=.6)
    # bootstrap clouds
    for Sb,c in ((S1b,'#d62728'),(S2b,'#2ca02c'),(S3b,'#1f77b4')):
        xs=[]; ys=[]
        for tr,pl in Sb:
            x,y=tp2xy(tr,pl); xs.append(x); ys.append(y)
        a.scatter(xs,ys,s=4,c=c,alpha=0.12,edgecolors='none')
    # P/T axes of individual events (data)
    for _,ev in q.iterrows():
        n,u=sdr_to_normal_slip(ev.s1,ev.d1,ev.r1)
        P=(u-n)/np.linalg.norm(u-n); Tx=(u+n)/np.linalg.norm(u+n)
        for V,mk,cc in ((P,'o','0.25'),(Tx,'s','0.55')):
            tr,pl=axis_trend_plunge(V); x,y=tp2xy(tr,pl)
            a.scatter([x],[y],s=26,marker=mk,facecolors='none',edgecolors=cc,lw=1.1)
    for S,c,lab in ((s1,'#d62728','σ1'),(s2,'#2ca02c','σ2'),(s3,'#1f77b4','σ3')):
        tr,pl=axis_trend_plunge(S); x,y=tp2xy(tr,pl)
        a.scatter([x],[y],s=150,c=c,edgecolors='k',lw=1.4,zorder=5)
        a.annotate(lab,(x,y),fontsize=13,fontweight='bold',ha='center',va='center',
                   color='white',zorder=6)
    a.text(0,1.08,'N',ha='center',fontsize=12,fontweight='bold')
    a.set_xlim(-1.15,1.15); a.set_ylim(-1.15,1.2); a.set_aspect('equal'); a.axis('off')
    a.set_title(f'Principal stress axes (equal-area, lower hemisphere)\n'
                f'n={len(q)} mechanisms · bootstrap clouds (open circ=P, open sq=T axes)',
                fontsize=10)
    # SHmax rose + R histogram
    b=ax[1]
    b.hist(Rb, bins=20, color='0.7', edgecolor='k')
    b.axvline(R, color='crimson', lw=2)
    b.set_xlabel('stress ratio R = (σ1−σ2)/(σ1−σ3)'); b.set_ylabel('bootstrap count')
    sh_mean, sh_sd = circ_axial(AZb)
    b.set_title(f'SHmax azimuth = {sh_mean:.0f}° ± {sh_sd:.0f}° (axial, 1σ)\n'
                f'R = {R:.2f}  (bootstrap {np.percentile(Rb,16):.2f}–{np.percentile(Rb,84):.2f})',
                fontsize=10)
    fig.tight_layout(); fig.savefig(fname, dpi=200); plt.close(fig)

def main(argv):
    def opt(flag, default, cast=float):
        return cast(argv[argv.index(flag)+1]) if flag in argv else default
    vrmin=opt('--vrmin',0.3); dcmin=opt('--dcmin',50); nmin=opt('--nmin',6,int)
    nusedmin=opt('--nusedmin',5,int); nboot=opt('--nboot',2000,int)
    alld,q = load_quality(vrmin,dcmin,nmin,nusedmin)
    print(f"catalog ok={len(alld)}  quality(VR>={vrmin},DC>={dcmin},nsta_used>={nusedmin})={len(q)}")
    if len(q) < 5:
        print(f"INSUFFICIENT: need >=5 (ideally >=12) mechanisms for a meaningful "
              f"inversion; have {len(q)}. Re-run after the batch fills in.")
        return 2
    sdr=[(r.s1,r.d1,r.r1) for _,r in q.iterrows()]
    nn=[]; uu=[]
    # plane-independent P/T axes (intrinsic to each moment tensor)
    Ptr=[]; Ppl=[]; Ttr=[]; Tpl=[]
    for st,dp,rk in sdr:
        n,u=sdr_to_normal_slip(st,dp,rk); nn.append(n); uu.append(u)
        P=(u-n)/np.linalg.norm(u-n); Tx=(u+n)/np.linalg.norm(u+n)
        tr,pl=axis_trend_plunge(P); Ptr.append(tr); Ppl.append(pl)
        tr,pl=axis_trend_plunge(Tx); Ttr.append(tr); Tpl.append(pl)
    Pm,Ps=circ_axial(Ptr); Tm,Ts=circ_axial(Ttr)
    pt=dict(P_trend=round(Pm,1), P_trend_sd=round(Ps,1), P_plunge_med=round(float(np.median(Ppl)),1),
            T_trend=round(Tm,1), T_trend_sd=round(Ts,1), T_plunge_med=round(float(np.median(Tpl)),1),
            SHmax_from_T=round((Tm+90)%180,1))   # SHmax perp to horizontal extension
    T,mis = invert_stress(nn,uu)
    s1,s2,s3,R = principal_axes(T)
    az,_ = shmax_azimuth(T)
    boot = bootstrap(sdr, nboot=nboot)
    t1,p1=axis_trend_plunge(s1); t2,p2=axis_trend_plunge(s2); t3,p3=axis_trend_plunge(s3)
    S1b,S2b,S3b,Rb,AZb = boot
    # confidence cones (deg) about each mean bootstrap axis
    cone1=axis_cone(S1b, mean_axis(S1b)); cone2=axis_cone(S2b, mean_axis(S2b))
    cone3=axis_cone(S3b, mean_axis(S3b))
    sh_mean, sh_sd = circ_axial(AZb)
    res=dict(n=len(q), vrmin=vrmin, dcmin=dcmin, nusedmin=nusedmin, mean_misfit_deg=round(float(mis),2),
             sigma1=dict(trend=round(t1,1),plunge=round(p1,1),cone68=round(cone1,1)),
             sigma2=dict(trend=round(t2,1),plunge=round(p2,1),cone68=round(cone2,1)),
             sigma3=dict(trend=round(t3,1),plunge=round(p3,1),cone68=round(cone3,1)),
             R=round(float(R),3), SHmax_az=round(float(az),1),
             SHmax_axial_mean=round(sh_mean,1), SHmax_axial_sd=round(sh_sd,1),
             R_CI68=[round(float(np.percentile(Rb,16)),3),round(float(np.percentile(Rb,84)),3)],
             n_SS=int(((np.abs(((q.r1+180)%360)-180)<30)|(np.abs(np.abs(((q.r1+180)%360)-180)-180)<30)).sum()),
             ptaxes=pt)
    os.makedirs(os.path.dirname(OUTJ),exist_ok=True)
    json.dump(res, open(OUTJ,'w'), indent=2)
    os.makedirs(os.path.dirname(OUTF),exist_ok=True)
    plot_stereonet(q, s1,s2,s3, boot, az, R, OUTF)
    print(json.dumps(res, indent=2)); print("wrote", OUTJ, "and", OUTF)
    return 0

if __name__=="__main__":
    if "--selftest" in sys.argv:
        sys.exit(0 if selftest() else 1)
    sys.exit(main(sys.argv))
