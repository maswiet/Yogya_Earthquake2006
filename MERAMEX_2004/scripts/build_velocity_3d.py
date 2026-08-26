#!/usr/bin/env python3
"""Inject the JSON payloads from gen_velocity_3d_data.py into velocity_3d.template.html
to produce a self-contained interactive artifact (no external JS libraries; the
template hand-rolls a small 3-D projection on a 2-D canvas).

Usage:
  python3 gen_velocity_3d_data.py --out ../figures/velocity_3d_data
  python3 build_velocity_3d.py --data ../figures/velocity_3d_data --out ../figures/velocity_3d.html
"""
import argparse, os

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--template", default=os.path.join(HERE, "velocity_3d.template.html"))
    ap.add_argument("--data", default=os.path.join(HERE, "..", "figures", "velocity_3d_data"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "figures", "velocity_3d.html"))
    a = ap.parse_args()

    with open(a.template, encoding="utf-8") as f:
        tpl = f.read()

    with open(os.path.join(a.data, "cells.json")) as f:
        cells = f.read()
    with open(os.path.join(a.data, "events.json")) as f:
        events = f.read()
    with open(os.path.join(a.data, "coast.json")) as f:
        coast = f.read()

    out = (tpl.replace("__CELLS_DATA__", cells)
              .replace("__EVENTS_DATA__", events)
              .replace("__COAST_DATA__", coast))

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(out)
    print(f"wrote {a.out} ({len(out)/1024:.1f} KB)")


if __name__ == "__main__":
    main()
