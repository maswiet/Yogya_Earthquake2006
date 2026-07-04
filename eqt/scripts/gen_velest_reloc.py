#!/usr/bin/env python3
"""Build NLLoc inputs that use the VELEST 1-D model + station corrections, to
relocate the full catalog. Parses the final VELEST P/S model and station
corrections from velest/yogya.OUT and writes:
  nll/nll_vel_v.in   (Vel2Grid, LAYER with VELEST Vp & Vs)
  nll/nll_time_v.in  (Grid2Time P and S)
  nll/nll_loc_v.in   (NLLoc with LOCDELAY station corrections, output yogya_v)
Reuses nll/obs/eqt.obs and nll/stations_gtsrce.txt.
"""
import os, re
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = f"{ROOT}/velest/yogya.OUT"
NLL = f"{ROOT}/nll"

lines = open(OUT).read().splitlines()
row = re.compile(r"^\s*([\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s*$")   # vel dv depth
def read_model(tag):
    st = max(i for i,l in enumerate(lines) if l.strip()==f"Velocity model   {tag}")
    blk=[]
    for l in lines[st+1:]:
        m=row.match(l)
        if not m: break
        blk.append((float(m.group(3)), float(m.group(1))))   # (depth, vel)
    return blk
P = dict(read_model(1)); S = dict(read_model(2))
depths = sorted(P)
print("VELEST model layers:", [(d,P[d],round(S[d],2)) for d in depths])

# station corrections (last 'Adjusted station corrections' block)
i = max(i for i,l in enumerate(lines) if "Adjusted station corrections:" in l)
cor={}
for l in lines[i+2:i+8]:
    t=l.split(); k=0
    while k+1<len(t):
        try: cor[t[k]]=float(t[k+1]); k+=3
        except (ValueError,IndexError): break
print(f"{len(cor)} station P-corrections")

HEAD="CONTROL 1 54321\nTRANS SIMPLE -7.92 110.44 0.0\n"

# --- Vel2Grid control (LAYER carries both Vp and Vs) ---
vg=[HEAD,"VGOUT ./model/layer_v","VGTYPE P","VGTYPE S",
    "VGGRID 111 111 51 -55.0 -55.0 -3.0 1.0 1.0 1.0 SLOW_LEN"]
# NLLoc/Vel2Grid: run twice (P then S); we write both VGTYPE and run selectively
lay=[]
for d in depths:
    lay.append(f"LAYER {d:6.1f}  {P[d]:.2f} 0.00  {S[d]:.2f} 0.00  2.70 0.00")
open(f"{NLL}/nll_vel_v.in","w").write(HEAD+"VGOUT ./model/layer_v\nVGTYPE P\n"+
    "VGGRID 111 111 51 -55.0 -55.0 -3.0 1.0 1.0 1.0 SLOW_LEN\n"+"\n".join(lay)+"\n")
open(f"{NLL}/nll_vel_v_S.in","w").write(HEAD+"VGOUT ./model/layer_v\nVGTYPE S\n"+
    "VGGRID 111 111 51 -55.0 -55.0 -3.0 1.0 1.0 1.0 SLOW_LEN\n"+"\n".join(lay)+"\n")

# --- Grid2Time controls (P and S) ---
gtsrce=open(f"{NLL}/stations_gtsrce.txt").read()
for ph in ("P","S"):
    open(f"{NLL}/nll_time_v_{ph}.in","w").write(
        HEAD+f"GTFILES ./model/layer_v ./time/layer_v {ph}\nGTMODE GRID3D ANGLES_NO\n"
        +gtsrce+"\nGT_PLFD 1.0e-3 0\n")

# --- NLLoc control with LOCDELAY station corrections ---
delay="\n".join(f"LOCDELAY {c:<6} P 1 {v:6.3f}" for c,v in cor.items())
loc=f"""{HEAD}LOCSIG "EQT Yogya 2006 - VELEST model reloc"
LOCFILES ./obs/eqt.obs NLLOC_OBS ./time/layer_v ./loc_v/yogya_v
LOCHYPOUT SAVE_NLLOC_ALL NLL_FORMAT_VER_2
LOCSEARCH OCT 12 12 8 0.01 20000 1000 0 1
LOCGRID 111 111 51 -55.0 -55.0 -3.0 1.0 1.0 1.0 PROB_DENSITY SAVE
LOCMETH EDT_OT_WT 9999.0 4 -1 -1 -1 -1 -1 1
LOCGAU 0.2 0.0
LOCGAU2 0.02 0.05 2.0
LOCPHASEID P P
LOCPHASEID S S
LOCQUAL2ERR 0.1 0.2 0.5 1.0 2.0
{delay}
LOCANGLES ANGLES_NO 5
"""
open(f"{NLL}/nll_loc_v.in","w").write(loc)
print("wrote nll_vel_v(.S).in, nll_time_v_{P,S}.in, nll_loc_v.in")
