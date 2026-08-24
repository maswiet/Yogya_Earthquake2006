#!/usr/bin/env python3
"""Checkerboard recovery against depth for the three data configurations.

One panel per checker size; correlation on the left axis, amplitude recovery on
the right. This is the figure that decides whether a LOTOS run is worth it and
what it should be run on.
"""
import argparse, collections, os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

LABEL = {"A_koulakov": "13,000 rays — Koulakov 2007 data volume",
         "B_land": "21,363 rays — this study, land only",
         "C_land_obs": "37,023 rays — land + ocean-bottom stations"}
COLOR = {"A_koulakov": "#888888", "B_land": "#1f77b4", "C_land_obs": "#d62728"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", required=True, help="tag cell depth n corrP ampP corrS ampS")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = [l.split() for l in open(a.table) if l.strip()]
    by = collections.defaultdict(list)
    for tag, cell, dep, n, cp, ap_, cs, as_ in rows:
        d0, d1 = (float(x) for x in dep.split("-"))
        by[(cell, tag)].append((0.5 * (d0 + d1), float(cp), float(ap_.rstrip("%")),
                                float(cs), float(as_.rstrip("%"))))

    cells = sorted({c for c, _ in by}, key=lambda x: -float(x))
    fig, axes = plt.subplots(1, len(cells), figsize=(6.2 * len(cells), 5.4), sharey=True)
    if len(cells) == 1:
        axes = [axes]
    for ax, cell in zip(axes, cells):
        for tag in ("A_koulakov", "B_land", "C_land_obs"):
            v = sorted(by.get((cell, tag), []))
            if not v:
                continue
            z = [x[0] for x in v]; c = [x[1] for x in v]
            ax.plot(c, z, "-o", ms=4, lw=1.8, color=COLOR[tag], label=LABEL[tag])
        ax.axvline(0.7, color="k", ls="--", lw=0.9)
        ax.text(0.705, 2, "conventionally resolved", rotation=90, va="bottom",
                fontsize=8, color="0.3")
        ax.set_xlim(0, 1); ax.set_ylim(65, 0)
        ax.set_xlabel("pattern correlation (P wave)")
        ax.set_title(f"{cell} km checkerboard", fontsize=12)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("depth (km)")
    axes[0].legend(loc="lower left", fontsize=8.5, framealpha=0.95)
    fig.suptitle("Checkerboard recovery vs depth — ±7 % anomalies, picking noise "
                 "0.1 s (P) / 0.2 s (S), origin times free", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(a.out, dpi=160, bbox_inches="tight")
    print("wrote", os.path.abspath(a.out))


if __name__ == "__main__":
    main()
