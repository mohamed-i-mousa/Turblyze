# Reproducing The Sphere Drag Validation

Commands assume the Turblyze executable is in `build.nosync/`.

## Mesh

Meshes are not shipped. Generate or copy the sphere Fluent mesh to:

```text
validation/turblyze/meshes/sphere.msh
```

## Laminar Sweep

```bash
cd build.nosync
CASE=../validation/turblyze/sphereLaminarCase
OUT=../validation/results/laminarCasesResults.csv
echo "Re,Cd,Cl,CdPressure,CdFriction" > "$OUT"

for Re in 1 2 5 10 20 50 100 200 500 1000 2000 5000; do
    mu=$(python3 -c "print(2.45/$Re)")
    sed -e "s|mu .*Pa.s.*|mu ${mu};|" \
        -e "s|filename .*\\.vtu;|filename /tmp/laminarRe${Re}.vtu;|" \
        "$CASE" > "/tmp/laminarRe${Re}.case"
    ./Turblyze "/tmp/laminarRe${Re}.case" > "/tmp/laminarRe${Re}.log" 2>&1
    f="/tmp/laminarRe${Re}_forces.txt"
    cdt=$(grep '^Cd' "$f" | awk '{print $4}')
    cdp=$(grep '^Cd' "$f" | awk '{print $2}')
    cdf=$(grep '^Cd' "$f" | awk '{print $3}')
    clt=$(grep '^Cl' "$f" | awk '{print $4}')
    echo "${Re},${cdt},${clt},${cdp},${cdf}" >> "$OUT"
done
```

## Subcritical SST Sweep

```bash
cd build.nosync
CASE=../validation/turblyze/sphereSSTCase
OUT=../validation/results/turbulentCasesResults.csv
echo "Re,Cd,Cl,CdPressure,CdFriction" > "$OUT"

for Re in 10000 50000 100000; do
    mu=$(python3 -c "print(2.45/$Re)")
    sed -e "s|mu .*Pa.s.*|mu ${mu};|" \
        -e "s|numIterations  *[0-9][0-9]*;|numIterations 1000;|" \
        -e "s|filename .*\\.vtu;|filename /tmp/turbRe${Re}.vtu;|" \
        "$CASE" > "/tmp/turbRe${Re}.case"
    ./Turblyze "/tmp/turbRe${Re}.case" > "/tmp/turbRe${Re}.log" 2>&1
    f="/tmp/turbRe${Re}_forces.txt"
    cdt=$(grep '^Cd' "$f" | awk '{print $4}')
    cdp=$(grep '^Cd' "$f" | awk '{print $2}')
    cdf=$(grep '^Cd' "$f" | awk '{print $3}')
    clt=$(grep '^Cl' "$f" | awk '{print $4}')
    echo "${Re},${cdt},${clt},${cdp},${cdf}" >> "$OUT"
done
```

The existing Re ~= 1.37e5 SST operating point is included in
`results/turbulentCasesResults.csv` as a derived total-Cd row.

## Plot

```bash
python3 validation/scripts/plotDragCurve.py
```
