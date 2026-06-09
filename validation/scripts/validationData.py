
import os

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))

def readCdData(filename):
    """
    Return (Re, Cd, Cl) lists from a drag CSV (laminar or turbulent).

    Columns: Re, Cd, Cl, CdPressure, CdFriction. Sorted by Re.
    """
    import csv
    rows = []
    with open(os.path.join(RESULTS, filename), newline="") as fh:
        for row in csv.DictReader(fh):
            cl = float(row["Cl"]) if row["Cl"] else float("nan")
            rows.append((float(row["Re"]), float(row["Cd"]), cl))
    rows.sort()
    reynolds = [r[0] for r in rows]
    cd = [r[1] for r in rows]
    cl = [r[2] for r in rows]
    return reynolds, cd, cl
