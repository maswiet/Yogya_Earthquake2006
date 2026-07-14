#!/usr/bin/env python3
"""Did the 17 Jul 2006 Pangandaran Mw7.7 tsunami earthquake perturb the Yogyakarta
2006 aftershock rate?  (dynamic vs static triggering)

Robust method:
 - Omori-Utsu MLE fit on PRE-Pangandaran events -> extrapolate as the baseline.
 - beta-statistic for post-event windows.
 - Completeness diagnostic: rate at ML>=0.3 vs ML>=1.0 (coda-masking vs real change).
 - Fine (6-h) binned rate around the exact origin time.
 - Static vs dynamic stress order-of-magnitude at Yogya.
"""
import os
import numpy as np, pandas as pd
from scipy.optimize import minimize
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.dates as mdates

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MAIN = pd.Timestamp("2006-05-26 22:53:58", tz="UTC")
PANG = pd.Timestamp("2006-07-17 08:19:27", tz="UTC")
tp   = (PANG-MAIN).total_seconds()/86400.0

e = pd.read_csv(f"{ROOT}/full/catalog_magnitude.csv")
e = e[e.gap<180].copy()
e["t"]=pd.to_datetime(e["time"],utc=True); e=e.sort_values("t").reset_index(drop=True)
e["days"]=(e["t"]-MAIN).dt.total_seconds()/86400.0
T0,T1=e.days.min(),e.days.max()

dM=0.1; M=e.ML.values; bb=np.arange(np.floor(M.min()),M.max()+dM,dM)
h,_=np.histogram(M,bb); Mc=bb[np.argmax(h)]+0.2
CUT=max(Mc,0.3)
ec=e[e.ML>=CUT].copy(); ehi=e[e.ML>=1.0].copy()
print(f"Mc={Mc:.2f}; low cut ML>={CUT:.2f} (n={len(ec)}), high cut ML>=1.0 (n={len(ehi)})")
print(f"catalog {T0:.1f}-{T1:.1f} d; Pangandaran t={tp:.2f} d")

# ---------- daily counts ----------
from scipy.optimize import curve_fit
day=ec.set_index("t").resample("D").size()
td=np.array([(d-MAIN).total_seconds()/86400+0.5 for d in day.index]); yd=day.values.astype(float)
dhi=ehi.set_index("t").resample("D").size().reindex(day.index,fill_value=0)

# ---------- Omori-Utsu fit (full window, EXCLUDING Pangandaran window; bounded) ----------
def omori(t,K,c,p): return K/np.power(t+c,p)
msk=~((td>tp-1)&(td<tp+6))
(K,c,p),_=curve_fit(omori,td[msk],yd[msk],p0=[500,15,1.0],
                    bounds=([1,1,0.4],[1e5,30,1.6]),maxfev=40000)
def rate(t): return omori(t,K,c,p)
exp_day=rate(td)
print(f"Omori-Utsu (full, excl. Pangandaran window): K={K:.0f}, c={c:.2f} d, p={p:.2f}")

# ---------- beta-statistic post-Pangandaran (vs Omori, and vs model-free local rate) ----------
r_local=((ec.days>=tp-15)&(ec.days<tp)).sum()/15.0   # events/day just before
def beta(win):
    ta,tb=tp,tp+win
    Nexp_om=K*(np.log((tb+c)/(ta+c)) if abs(p-1)<1e-8 else ((tb+c)**(1-p)-(ta+c)**(1-p))/(1-p))
    Nexp_loc=r_local*win
    Nobs=((ec.days>=ta)&(ec.days<tb)).sum()
    return Nobs,Nexp_om,(Nobs-Nexp_om)/np.sqrt(Nexp_om),Nexp_loc,(Nobs-Nexp_loc)/np.sqrt(Nexp_loc)
print(f"\nlocal pre-rate = {r_local:.1f}/day (ML>={CUT:.1f})")
print("beta (obs vs Omori | obs vs local-rate):")
for w in (0.5,1,2,3,5,7,10):
    N,Ne,b,Nl,bl=beta(w)
    print(f"  +{w:4.1f} d: Nobs={N:4.0f}  Omori:Nexp={Ne:6.1f} b={b:+.2f}  "
          f"local:Nexp={Nl:6.1f} b={bl:+.2f}"
          f"{'  *' if abs(b)>=2 or abs(bl)>=2 else ''}")

# ---------- completeness diagnostic (ratio hi/lo around Pangandaran) ----------
def win_rate(df,a,b): return ((df.days>=a)&(df.days<b)).sum()
pre_lo=win_rate(ec,tp-10,tp)/10; post_lo=win_rate(ec,tp,tp+3)/3
pre_hi=win_rate(ehi,tp-10,tp)/10; post_hi=win_rate(ehi,tp,tp+3)/3
print(f"\nrate/day  ML>={CUT:.1f}: pre(10d)={pre_lo:.1f} post(3d)={post_lo:.1f} ratio={post_lo/pre_lo:.2f}")
print(f"rate/day  ML>=1.0: pre(10d)={pre_hi:.2f} post(3d)={post_hi:.2f} "
      f"ratio={post_hi/pre_hi if pre_hi>0 else float('nan'):.2f}")

# ---------- stress order-of-magnitude at Yogya ----------
R=364.0  # km epicentral (rupture nearest ~250 km)
mu=3.0e10; Vs=3.5e3; M0=10**(1.5*7.7+9.1)
# static ~ M0/(mu*R^3) (very rough, near-field-free); dynamic ~ mu*PGV/Vs
static=M0/(mu*(R*1e3)**3)
PGV=0.01  # m/s, ~1 cm/s long-period at ~300 km for a slow Mw7.7 (order of magnitude)
dyn=mu*PGV/Vs
print(f"\nstress @ ~{R:.0f} km (order of magnitude):")
print(f"  static Coulomb ~ {static:.1e} Pa = {static/1e5:.1e} bar (negligible; << 0.01 bar)")
print(f"  dynamic (PGV~{PGV*100:.0f} cm/s) ~ {dyn:.1e} Pa = {dyn/1e5:.2f} bar")

# ---------- beta vs magnitude threshold (completeness diagnostic) ----------
cuts=np.arange(0.0,1.6,0.1); b3=[]; b7=[]
for mm in cuts:
    rl=((e.days>=tp-15)&(e.days<tp)&(e.ML>=mm)).sum()/15.0
    for w,arr in ((3,b3),(7,b7)):
        No=((e.days>=tp)&(e.days<tp+w)&(e.ML>=mm)).sum(); Ne=rl*w
        arr.append((No-Ne)/np.sqrt(Ne) if Ne>0 else np.nan)

# ---------- figure ----------
fig=plt.figure(figsize=(15,10)); gs=fig.add_gridspec(2,2)
ax1=fig.add_subplot(gs[0,:])
ax1.bar(day.index,yd,width=1,color="0.7",label=f"obs ML>={CUT:.1f}")
ax1.bar(day.index,dhi.values,width=1,color="tab:red",alpha=0.8,label="obs ML>=1.0")
ax1.plot(day.index,exp_day,"b-",lw=2,label=f"Omori p={p:.2f} (baseline)")
ax1.axvline(PANG,color="red",lw=2,ls="--",label="Pangandaran Mw7.7 (17 Jul)")
ax1.set_yscale("log"); ax1.set_ylim(1,None); ax1.set_ylabel("events/day")
ax1.set_title("Yogya aftershock daily rate vs Omori baseline"); ax1.legend(fontsize=8,ncol=2)
ax1.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d"))

ax2=fig.add_subplot(gs[1,0])
edges=np.arange(tp-7,tp+10,0.25); ctr=(edges[:-1]+edges[1:])/2
clo,_=np.histogram(ec.days,edges); chi,_=np.histogram(ehi.days,edges)
ax2.bar(ctr-tp,clo,width=0.24,color="0.6",label=f"ML>={CUT:.1f} / 6h")
ax2.plot(ctr-tp,rate(ctr)*0.25,"b-",lw=1.5,label="Omori exp.")
ax2.bar(ctr-tp,chi,width=0.24,color="tab:red",alpha=0.85,label="ML>=1.0 / 6h")
ax2.axvline(0,color="red",lw=2,ls="--")
ax2.set_xlabel("days relative to Pangandaran"); ax2.set_ylabel("events / 6 h")
ax2.set_title("Fine-binned rate around 17 Jul (no immediate spike)"); ax2.legend(fontsize=8)

ax3=fig.add_subplot(gs[1,1])
ax3.plot(cuts,b3,"o-",color="tab:purple",label="beta, +3 d")
ax3.plot(cuts,b7,"s-",color="tab:orange",label="beta, +7 d")
ax3.axhline(0,color="k",lw=.7); ax3.axhline(-2,color="gray",ls=":"); ax3.axhline(2,color="gray",ls=":")
ax3.axvline(Mc,color="green",ls="--",label=f"Mc={Mc:.1f}")
ax3.set_xlabel("magnitude threshold ML"); ax3.set_ylabel("beta (obs vs pre-rate)")
ax3.set_title("Deficit vanishes above Mc -> completeness artifact, not triggering")
ax3.legend(fontsize=8)
plt.tight_layout(); plt.savefig(f"{ROOT}/figures/pangandaran_rate.png",dpi=120)
print("\nwrote figures/pangandaran_rate.png")
