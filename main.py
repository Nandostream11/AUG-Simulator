import argparse
from Modeling.glider_model import Vertical_Motion
from Modeling.dynamics import Dynamics
from Parameters.slocum import SLOCUM_PARAMS
import numpy as np


def main(args):
    Z = Vertical_Motion(args)
    Z.set_desired_trajectory()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="An Autonomous Underwater Glider Simulator."
    )
    parser.add_argument(
        "-i", "--info", help="give full information in each cycle", action="store_true"
    )
    parser.add_argument(
        "-c",
        "--cycle",
        help="number of desired cycles in sawtooth trajectory",
        default=4,
        type=int,
    )
    parser.add_argument(
        "-g",
        "--glider",
        help="desired glider mode.",
        default=("slocum"),
    )
    parser.add_argument(
        "-a",
        "--angle",
        help="desired glider angle",
        default=SLOCUM_PARAMS.VARIABLES.GLIDE_ANGLE,
        type=int,
    )
    parser.add_argument(
        "-s",
        "--speed",
        help="desired glider speed",
        default=SLOCUM_PARAMS.VARIABLES.SPEED,
        type=float,
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="variables to be plotted [all, x, z, omega, v, rp1, rp3, mb, theta]",
        default="all",
        nargs="*",
    )

    args = parser.parse_args()

    main(args)
