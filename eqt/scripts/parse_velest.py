#!/usr/bin/env python3
"""Parse the final VELEST 1-D model + station corrections from yogya.OUT and plot
against the starting model. Writes velest/velest_result.png and prints the model.
"""
import os, re
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = f"{ROOT}/velest/yogya.OUT"

# ---- parse the last 'Velocity model 1' (P) and 'Velocity model 2' (S) blocks ----
lines_all = open(OUT).read().splitlines()
row = re.compile(r"^\s*([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")   # vel  dv  depth
def read_model(tag):
    starts = [i for i,l in enumerate(lines_all) if l.strip()==f"Velocity model   {tag}"]
    blk=[]
    for l in lines_all[starts[-1]+1:]:
        m=row.match(l)
        if not m: break
        vel=float(m.group(1)); dep=float(m.group(3))
        blk.append((dep, vel))
    return blk
pblk = read_model(1); sblk = read_model(2)
print("FINAL VELEST 1-D MODEL")
print(" depth(km)   Vp     Vs    Vp/Vs")
sd = dict(sblk)
for (z, vp) in pblk:
    vs = sd.get(z, np.nan)
    print(f"  {z:6.1f}   {vp:5.2f}  {vs:5.2f}   {vp/vs:5.3f}" if vs==vs else f"  {z:6.1f}   {vp:5.2f}")

# ---- station corrections (last 'Adjusted station corrections' block) ----
lines = open(OUT).read().splitlines()
idx = max(i for i,l in enumerate(lines) if "Adjusted station corrections:" in l)
cor = {}
for l in lines[idx+2:idx+8]:
    toks = l.split()
    k = 0
    while k+1 < len(toks):
        try:
            cor[toks[k]] = float(toks[k+1]); k += 3
        except (ValueError, IndexError):
            break
print(f"\nstation P-corrections: {len(cor)} stations, "
      f"range {min(cor.values()):+.2f}..{max(cor.values()):+.2f} s")

# ---- starting model (from gen_velest) ----
START = [(-2.0,4.00),(0.0,4.50),(3.0,5.20),(6.0,5.70),(9.0,6.00),(12.0,6.20),
         (16.0,6.40),(20.0,6.70),(28.0,7.20),(40.0,8.00)]

def step(mod):
    z=[]; v=[]
    for i,(d,vv) in enumerate(mod):
        z.append(d); v.append(vv)
    return np.array(z), np.array(v)

fig, ax = plt.subplots(1, 2, figsize=(13, 6))
zs, vps = step(START)
ax[0].step(vps, zs, where="post", color="gray", ls="--", label="start Vp")
zp = np.array([z for z,_ in pblk]); vp = np.array([v for _,v in pblk])
zssm = np.array([z for z,_ in sblk]); vs = np.array([v for _,v in sblk])
ax[0].step(vp, zp, where="post", color="tab:blue", lw=2, label="VELEST Vp")
ax[0].step(vs, zssm, where="post", color="tab:red", lw=2, label="VELEST Vs")
ax[0].set_ylim(35, -3); ax[0].set_xlabel("velocity (km/s)"); ax[0].set_ylabel("depth (km)")
ax[0].set_title("VELEST minimum-1D model"); ax[0].legend(); ax[0].grid(alpha=0.3)

names = list(cor.keys()); vals = [cor[n] for n in names]
order = np.argsort(vals)
ax[1].barh([names[i] for i in order], [vals[i] for i in order], color="teal")
ax[1].set_xlabel("P station correction (s)"); ax[1].set_title("Station corrections (relative)")
ax[1].axvline(0, color="k", lw=0.8); ax[1].grid(alpha=0.3, axis="x")
plt.tight_layout(); plt.savefig(f"{ROOT}/velest/velest_result.png", dpi=110)
print("wrote velest/velest_result.png")
