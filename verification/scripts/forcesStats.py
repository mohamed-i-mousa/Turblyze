#!/usr/bin/env python3
"""
Mean/std of Cd and Cl over the limit cycle from a force-history CSV.

The <case>_forces_history.csv input (one row per SIMPLE iteration) is a
retained verification artifact; the current solver writes only a converged
<case>_forces.txt summary, not a per-iteration history. The sphere wake at
Re=1.37e5 is unsteady, so the steady solve settles into a statistically-steady
limit cycle rather than a fixed point; the reported coefficient is the mean
over that cycle (after the initial transient).

Usage:
    python3 forcesStats.py ../results/sphereForcesHistory.csv
    python3 forcesStats.py ../results/sphereForcesHistory.csv --start 1000
"""

import argparse
import csv


def readHistory(path):
    rows = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append((int(row["iteration"]), float(row["Cd"]),
                         float(row["Cl"])))
    return rows


def meanStd(values):
    n = len(values)
    mean = sum(values) / n
    var = sum((v - mean) ** 2 for v in values) / n
    return mean, var ** 0.5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("history", help="*_forces_history.csv")
    parser.add_argument("--start", type=int, default=None,
                        help="first iteration to average (default: second half)")
    args = parser.parse_args()

    rows = readHistory(args.history)
    last = rows[-1][0]
    start = args.start if args.start is not None else last // 2
    window = [(cd, cl) for (it, cd, cl) in rows if it >= start]

    cdMean, cdStd = meanStd([cd for cd, _ in window])
    clMean, clStd = meanStd([cl for _, cl in window])

    print(f"file       : {args.history}")
    print(f"window     : iters {start}-{last}  (n={len(window)})")
    print(f"Cd  mean   : {cdMean:.4f}  std {cdStd:.4f}")
    print(f"Cl  mean   : {clMean:.4f}  std {clStd:.4f}")


if __name__ == "__main__":
    main()
