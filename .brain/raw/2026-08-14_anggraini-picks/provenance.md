# Provenance — Anggraini (2013) manual pick data

- Source: user-provided (the group's own dissertation data), 2026-08-14.
- Files (immutable copies):
  - `phase_300.dat` — 588 events, 3776 P + 3414 S manual picks. Per-event header
    `# y mo d h mi s lat lon dep ML ... evid`, then phase lines
    `STA  travel_time_s  weight  P|S`. Covers 3-7 June 2006 (full 24h origin times).
  - `station.dat` — manual station codes + lat/lon (YOG,TRI,PRA,DES,KRI,WON,PEL,
    RAT,NGL,WAN,BOG,BUM,PAL).
  - `stat_ft.dat` — 10 stations, lon/lat + P station correction (s) + numbered code.
- One bad row (ML 10.89 / negative timestamp) filtered in analysis.
- This is the PICK-level ground truth referenced in Diambama et al. 2019
  (3769 P + 3407 S, 588 events). Origin times here (h:m:s) are what the
  Bantul2006 xlsx catalogue lacked.
