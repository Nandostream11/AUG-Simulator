import argparse
from Modeling2d.glider_model_2D import Vertical_Motion
from Modeling3d.glider_model_3D import ThreeD_Motion
from Parameters.slocum import SLOCUM_PARAMS


def main(args):
    if args.mode == "2D":
        Z = Vertical_Motion(args)
    elif args.mode == "3D":
        Z = ThreeD_Motion(args)
    Z.set_desired_trajectory()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="An Autonomous Underwater Glider Simulator."
    )
    parser.add_argument(
        "-i", "--info", help="give full information in each cycle", action="store_true"
    )
    parser.add_argument("-m", "--mode", help="set mode as 2D or 3D", default="2D")
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
        help="desired glider model ['slocum']",
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
        "-pid",
        "--pid",
        help="enable or disable PID pitch control",
        default="disable",
    )
    parser.add_argument(
        "-p",
        "--plot",
        help="variables to be plotted [3D, all, x, y, z, omega1, omega2, omega3, vel, v1, v2, v3, rp1, rp2, rp3, mb, phi, theta, psi]",
        default="all",
        nargs="*",
    )

    args = parser.parse_args()

    main(args)
