#!/usr/bin/env python3
"""Build STATIC BayesISOLA inputs for Yogya 2006 MT inversion:
  eqt/mt/stations.xml  - ObsPy inventory, L4-3D response for all recorders (YK:TF07..TF19, HHZ/N/E)
  eqt/mt/crustal.dat   - refined VELEST 1-D model (depth Vp Vs Rho Qp Qs)
Response identical for all (L4-3D geophone, 1 Hz, h=0.707, 1.7e8 counts/(m/s), M/S).
Coordinates are supplied per-event via network.stn (BayesISOLA reads geometry there, response from XML).
"""
import os, json, numpy as np
from obspy import UTCDateTime
from obspy.core.inventory import Inventory, Network, Station, Channel, Site
from obspy.core.inventory.response import Response

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
MT = f"{ROOT}/mt"; os.makedirs(MT, exist_ok=True)
per = json.load(open(f"{ROOT}/config/stations_periods.json"))

# ---- L4-3D PAZ (velocity geophone): 2 zeros at 0; poles from f0=1 Hz, h=0.707 ----
f0, h = 1.0, 0.707
w0 = 2*np.pi*f0
poles = [complex(-h*w0,  w0*np.sqrt(1-h*h)), complex(-h*w0, -w0*np.sqrt(1-h*h))]
zeros = [0+0j, 0+0j]
SENS = 1.7e8          # counts/(m/s)
f_ref = 5.0           # Hz (flat band reference)
resp = Response.from_paz(zeros, poles, stage_gain=SENS,
                         stage_gain_frequency=f_ref,
                         input_units="M/S", output_units="COUNTS",
                         normalization_frequency=f_ref,
                         pz_transfer_function_type="LAPLACE (RADIANS/SECOND)")

chans_def = [("HHZ",0.0,-90.0),("HHN",0.0,0.0),("HHE",90.0,0.0)]
stations=[]
for code,info in per.items():
    s0 = info["sites"][0]
    chs=[]
    for ch,az,dip in chans_def:
        chs.append(Channel(code=ch, location_code="", latitude=s0["lat"], longitude=s0["lon"],
                           elevation=0.0, depth=0.0, azimuth=az, dip=dip,
                           sample_rate=100.0, response=resp))
    stations.append(Station(code=code, latitude=s0["lat"], longitude=s0["lon"], elevation=0.0,
                            site=Site(name=code), channels=chs))
inv = Inventory([Network(code="YK", stations=stations, description="GFZ XN/YK Yogya 2006 temp array (L4-3D)")],
                source="Yogya2006-MT")
inv.write(f"{MT}/stations.xml", format="STATIONXML")
print(f"wrote {MT}/stations.xml : {len(stations)} stations, poles={poles}")

# ---- crustal model: refined VELEST, density via Brocher(2005), Q crustal ----
LAYERS=[(0.0,2.90),(0.7,4.30),(2.0,4.65),(4.0,5.49),(10.0,6.30),
        (13.0,6.39),(16.0,6.55),(22.0,6.80),(26.0,7.20),(40.0,8.00)]
VPVS=1.735
def brocher_rho(vp): # g/cm^3
    return (1.6612*vp -0.4721*vp**2 +0.0671*vp**3 -0.0043*vp**4 +0.000106*vp**5)
with open(f"{MT}/crustal.dat","w") as f:
    f.write("Crustal model                refined VELEST (Yogya 2006)\n")
    f.write("number of layers\n   {0:d}\n".format(len(LAYERS)))
    f.write("Parameters of the layers\ndepth of layer top(km)   Vp(km/s)    Vs(km/s)    Rho(g/cm**3)    Qp     Qs\n")
    for d,vp in LAYERS:
        vs=vp/VPVS; rho=brocher_rho(vp); qp=300.0 if vp<6.0 else 500.0; qs=qp/2
        f.write(f"  {d:6.2f}   {vp:7.3f}   {vs:7.3f}   {rho:7.3f}   {qp:6.0f} {qs:6.0f}\n")
print(f"wrote {MT}/crustal.dat : {len(LAYERS)} layers, Vp/Vs={VPVS}")
