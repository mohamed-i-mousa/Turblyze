<!--
SPDX-FileCopyrightText: 2025-2026 Mohamed Mousa
SPDX-License-Identifier: Apache-2.0
-->

# Reproducing The OpenFOAM Verification

Commands assume the Turblyze executable is in `build.nosync/`.

## Mesh

Meshes are not shipped. Generate or copy the sphere Fluent mesh to:

```text
verification/turblyze/meshes/sphere.msh
```

Convert the same mesh into the OpenFOAM case from `verification/openFoam/`:

```bash
fluent3DMeshToFoam ../turblyze/meshes/sphere.msh
```

## Turblyze

```bash
cd build.nosync
./Turblyze ../verification/turblyze/sphereSSTCase

cp ../verification/results/sphere_forces.txt \
   ../verification/results/turblyzeForcesFull.txt
```

The current build writes only the converged `sphere_forces.txt` summary shown
above. The limit-cycle `Cd`/`Cl` mean and std reported in `README.md` were
computed by `forcesStats.py` over the per-iteration force history from the
original verification run, retained as the committed artifact
`verification/results/sphereForcesHistory.csv` (the current solver does not
regenerate it):

```bash
cd ..
python3 verification/scripts/forcesStats.py \
    verification/results/sphereForcesHistory.csv
```

## OpenFOAM

```bash
cd verification/openFoam
foamRun
cp postProcessing/forceCoeffs1/*/coefficient.dat \
   ../results/openFoamForceCoeffsFull.dat
```

`results/openFoamSummary.csv` stores the compact OpenFOAM summary used in the
README table.
