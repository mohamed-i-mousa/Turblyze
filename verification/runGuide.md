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
cp ../verification/results/sphere_forces_history.csv \
   ../verification/results/sphereForcesHistory.csv
```

Compute the Turblyze mean over the trailing half:

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
