<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

# Verification: Turblyze vs OpenFOAM

This folder verifies Turblyze's k-omega SST implementation against OpenFOAM on the same sphere case: same mesh, boundary conditions, model constants, reference velocity, and force-coefficient definition.

## Case

| Quantity | Value |
|---|---:|
| Sphere diameter | 0.100 m |
| Free-stream speed | 20.0 m/s |
| Density | 1.225 kg/m^3 |
| Dynamic viscosity | 1.7894e-5 Pa.s |
| Reynolds number | 1.37e5 |
| Turbulence model | k-omega SST |
| Momentum convection | Turblyze `SecondOrderUpwind`; OpenFOAM `linearUpwind` |

## Result

| Quantity | Turblyze | OpenFOAM | Difference |
|---|---:|---:|---:|
| Cd mean | 0.4648 | 0.5071 | 8.4% |
| Cd std | 0.0145 | 0.0002 | — |
| Cl mean | -0.0725 | -0.1417 | qualitative |

The drag agreement is within acceptable tolerance for this steady-RANS sphere case.

> **Force-statistics provenance**: the limit-cycle `Cd`/`Cl` mean and std were
> computed by `scripts/forcesStats.py` over the per-iteration force history from
> the original verification run, retained as the committed artifact
> `results/sphereForcesHistory.csv`. The current solver writes only a converged
> `<case>_forces.txt` summary, not a per-iteration history — see `runGuide.md`.

## Files

```text
verification/
  turblyze/sphereSSTCase
  openFoam/
  results/
    sphereForcesHistory.csv
    openFoamForceCoeffsFull.dat
    turblyzeForcesFull.txt
    openFoamSummary.csv
  scripts/forcesStats.py
```

See `runGuide.md` to reproduce the comparison.
