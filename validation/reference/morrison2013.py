#!/usr/bin/env python3
"""
Sphere drag coefficient from the Morrison (2013) correlation.

Single full-range fit Cd(Re) for a smooth sphere, valid for Re up to ~1e6,
including the drag-crisis region. Used as the experimental reference target
for the Turblyze sphere validation case.

Reference:
    Morrison, F. A. (2013). An Introduction to Fluid Mechanics.
    Cambridge University Press. (Sphere drag correlation, Fig. 8.13 / Eq.)

Usage:
    python3 morrison2013.py --Re 1.37e5
    python3 morrison2013.py --U 20 --D 0.1 --nu 1.461e-5
"""

import argparse

# Turblyze sphere case defaults
U_DEFAULT = 20.0          # freestream speed [m/s]
D_DEFAULT = 0.10          # sphere diameter [m]
NU_DEFAULT = 1.461e-5     # kinematic viscosity [m^2/s]


def morrisonCd(reynolds):
    """Drag coefficient of a smooth sphere at Reynolds number `reynolds`."""
    return (
        24.0 / reynolds
        + (2.6 * (reynolds / 5.0)) / (1.0 + (reynolds / 5.0) ** 1.52)
        + (0.411 * (reynolds / 2.63e5) ** -7.94)
        / (1.0 + (reynolds / 2.63e5) ** -8.00)
        + (0.25 * (reynolds / 1.0e6)) / (1.0 + (reynolds / 1.0e6))
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--Re", type=float, default=None,
                        help="Reynolds number (overrides U/D/nu)")
    parser.add_argument("--U", type=float, default=U_DEFAULT)
    parser.add_argument("--D", type=float, default=D_DEFAULT)
    parser.add_argument("--nu", type=float, default=NU_DEFAULT)
    args = parser.parse_args()

    reynolds = args.Re if args.Re is not None else args.U * args.D / args.nu
    cd = morrisonCd(reynolds)

    print(f"Re = {reynolds:.4g}")
    print(f"Cd (Morrison 2013) = {cd:.4f}")


if __name__ == "__main__":
    main()
