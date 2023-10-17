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
        "-i", "--info", help="Give full information in each cycle.", action="store_true"
    )
    parser.add_argument(
        "-c",
        "--cycle",
        help="Number of desired cycles in sawtooth trajectory.",
        default=4,
        type=int,
    )
    parser.add_argument(
        "-g",
        "--glider",
        help="Desired Glider Model.",
        default=("slocum"),
    )
    parser.add_argument(
        "-a",
        "--angle",
        help="Desired Glider Angle.",
        default=SLOCUM_PARAMS.VARIABLES.GLIDE_ANGLE,
        type=int,
    )
    parser.add_argument(
        "-s",
        "--speed",
        help="Desired Glider Speed.",
        default=SLOCUM_PARAMS.VARIABLES.SPEED,
        type=float,
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="Variable to be plotted. Values accepted are [all, x, z, omega, v, rp1, rp3, mb, theta].",
        default="all",
        nargs="*",
    )

    args = parser.parse_args()

    main(args)
