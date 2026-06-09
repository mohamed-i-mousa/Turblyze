#!/usr/bin/env python3
"""
Sphere drag curve: Turblyze vs Morrison (2013).

Two behaviours on one validation figure:

  * Laminar low-Re (turbulence OFF).
  * Subcritical k-omega SST points.

Writes figures/dragCurve.png.
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import validationData

sys.path.insert(0, os.path.join(validationData.HERE, "..", "reference"))
from morrison2013 import morrisonCd  # noqa: E402

FIG = os.path.join(validationData.HERE, "..", "figures", "dragCurve.png")


def main():
    laminarRe, laminarCd, _ = validationData.readCdData(
        "laminarCasesResults.csv"
    )
    turbulentRe, turbulentCd, _ = validationData.readCdData(
        "turbulentCasesResults.csv"
    )

    reynolds = np.logspace(np.log10(0.6), 6, 700)
    cd = [morrisonCd(r) for r in reynolds]

    fig, ax = plt.subplots(figsize=(8.5, 5.5))

    ax.loglog(reynolds, cd, color="black", lw=1.6,
              label="Morrison 2013 (Experiment)")

    # Laminar sweep (turbulence off) — validation.
    ax.loglog(laminarRe, laminarCd, "o", color="#1f77b4", ms=8, mec="black",
              mew=0.6, label="Turblyze — Laminar", zorder=5)

    # Subcritical fully-turbulent k-omega SST benchmark points.
    ax.loglog(turbulentRe, turbulentCd, "D", color="#d62728", ms=8, mec="black",
              mew=0.6, label="Turblyze — k-$\\omega$ SST", zorder=5)

    ax.set_xlabel("Reynolds number  Re = U·D/ν")
    ax.set_ylabel(r"$C_d$")
    ax.set_title("Sphere Drag Coefficient: Turblyze vs Experiment")
    ax.set_xlim(0.6, 1e6)
    ax.set_ylim(0.05, 100)
    ax.legend(fontsize=8.5, loc="lower left")
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    fig.savefig(FIG, dpi=600)
    print("wrote", os.path.normpath(FIG))


if __name__ == "__main__":
    main()
