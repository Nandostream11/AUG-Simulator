import numpy as np
import math
from scipy.integrate import solve_ivp
import utils
from Modeling.dynamics import Dynamics


class Vertical_Motion:
    def __init__(self, args):
        self.args = args
        self.cycles = self.args.cycle
        self.glider_name = self.args.glider
        self.info = self.args.info
        self.plots = self.args.plot

        self.initialization()

        self.solver_array = []
        self.total_time = []
        self.derivative = []

    def initialization(self):
        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        if self.glider_name == ("slocum"):
            from Parameters.slocum import SLOCUM_PARAMS as P
        else:
            print("Invalid glider model")
            raise ImportError

        self.mass_params = P.GLIDER_CONFIG
        self.hydro_params = P.HYDRODYNAMICS
        self.vars = P.VARIABLES

        self.mh = self.mass_params.HULL_MASS
        self.mw = self.mass_params.FIXED_POINT_MASS
        self.mb = self.mass_params.BALLAST_MASS
        self.mm = self.mass_params.INT_MOVABLE_MASS

        self.mt = self.mh + self.mw + self.mb + self.mm

        self.ms = self.mh + self.mw + self.mb
        self.mt = self.ms + self.mm

        self.m = self.mass_params.FLUID_DISP_MASS
        self.m0 = self.mt - self.m

        self.Mf = np.diag(
            [self.mass_params.MF1, self.mass_params.MF2, self.mass_params.MF3]
        )
        self.Jf = np.diag(
            [self.mass_params.J1, self.mass_params.J2, self.mass_params.J3]
        )

        self.M = self.mh * self.I3 + self.Mf
        self.J = self.Jf  # J = Jf + Jh

        self.KL = self.hydro_params.KL
        self.KL0 = self.hydro_params.KL0
        self.KD = self.hydro_params.KD
        self.KD0 = self.hydro_params.KD0
        self.KM = self.hydro_params.KM
        self.KM0 = self.hydro_params.KM0
        self.KOmega1 = self.hydro_params.KOmega1
        self.KOmega2 = self.hydro_params.KOmega2

        self.rp3 = self.vars.rp3
        self.rb1 = self.vars.rb1
        self.rb3 = self.vars.rb3

        # [self.rw1, self.rw2, self.rw3] = [0.0, 0.0, 0.0]

        self.glide_angle_deg = self.args.angle
        self.V_d = self.args.speed
        self.ballast_rate = self.vars.BALLAST_RATE

        self.set_first_run_params()

    def set_first_run_params(self):
        self.phi = self.vars.PHI
        self.theta0 = math.radians(self.vars.THETA)
        self.psi = self.vars.PSI

    def set_desired_trajectory(self):
        self.E_i_d = np.array(
            [
                math.radians(math.pow(-1, k + 1) * self.glide_angle_deg)
                for k in range(self.cycles)
            ]
        )

        self.lim1 = math.degrees(
            math.atan(
                2
                * (self.KD / self.KL)
                * (
                    self.KL0 / self.KL
                    + math.sqrt(math.pow(self.KL0 / self.KL, 2) + self.KD0 / self.KD)
                )
            )
        )

        self.lim2 = math.degrees(
            math.atan(
                2
                * (self.KD / self.KL)
                * (
                    self.KL0 / self.KL
                    - math.sqrt(math.pow(self.KL0 / self.KL, 2) + self.KD0 / self.KD)
                )
            )
        )

        l = len(self.E_i_d)
        for i in range(l):
            self.e_i_d = self.E_i_d[i]

            print(
                "\nIteration {} | Desired glide angle in deg = {}".format(
                    i, math.degrees(self.e_i_d)
                )
            )

            if (self.e_i_d) > 0:
                self.glider_direction = "U"
                self.ballast_rate = -abs(self.ballast_rate)
                print("Glider moving in upward direction")

            elif (self.e_i_d) < 0:
                self.glider_direction = "D"
                self.ballast_rate = abs(self.ballast_rate)
                print("Glider moving in downward direction")

            self.alpha_d = (
                (1 / 2)
                * (self.KL / self.KD)
                * math.tan(self.e_i_d)
                * (
                    -1
                    + math.sqrt(
                        1
                        - 4
                        * (self.KD / math.pow(self.KL, 2))
                        * (1 / math.tan(self.e_i_d))
                        * (self.KD0 * (1 / math.tan(self.e_i_d)) + self.KL0)
                    )
                )
            )

            self.mb_d = (self.m - self.mh - self.mm) + (1 / self.g) * (
                -math.sin(self.e_i_d) * (self.KD0 + self.KD * math.pow(self.alpha_d, 2))
                + math.cos(self.e_i_d) * (self.KL0 + self.KL * self.alpha_d)
            ) * math.pow(self.V_d, 2)

            self.m0_d = self.mb_d + self.mh + self.mm - self.m

            self.theta_d = self.e_i_d + self.alpha_d

            self.v1_d = self.V_d * math.cos(self.alpha_d)
            self.v3_d = self.V_d * math.sin(self.alpha_d)

            self.Pp1_d = self.mm * self.v1_d
            self.Pp3_d = self.mm * self.v3_d

            self.rp1_d = -self.rp3 * math.tan(self.theta_d) + (
                1 / (self.mm * self.g * math.cos(self.theta_d))
            ) * (
                (self.Mf[2, 2] - self.Mf[0, 0]) * self.v1_d * self.v3_d
                + (self.KM0 + self.KM * self.alpha_d) * math.pow(self.V_d, 2)
            )

            if self.info == True:
                print(
                    "Desired angle of attack in deg = {}".format(
                        math.degrees(self.alpha_d)
                    )
                )
                print("Desired ballast mass in kg = {}".format(self.mb_d))
                print(
                    "Desired position of internal movable mass in cm = {}".format(
                        self.rp1_d * 100
                    )
                )

            self.save_json()

            # Initial conditions at every peak of the sawtooth trajectory

            if i == 0:
                self.z_in = np.array(
                    [
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [self.v1_d, 0.0, self.v3_d],
                        [self.rp1_d, 0.0, self.rp3],
                        [self.rb1, 0.0, self.rb3],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [self.mb_d, 0, 0],
                        [self.theta0, 0, 0],
                    ]
                ).flatten()

            else:
                p = len(self.solver_array)
                self.z_in = np.array(
                    [
                        [self.solver_array[-1][0:3]],
                        [self.solver_array[-1][3:6]],
                        [self.solver_array[-1][6:9]],
                        [self.solver_array[-1][9:12]],
                        [self.solver_array[-1][12:15]],
                        [self.solver_array[-1][15:18]],
                        [self.solver_array[-1][18:21]],
                        [self.solver_array[-1][21:24]],
                        [self.solver_array[-1][24:27]],
                    ]
                ).flatten()

            self.i = i

            sol = self.solve_ode(self.z_in)

            if i == 0:
                self.solver_array = np.array(sol.y).T
                self.total_time = np.array(sol.t)
            else:
                self.solver_array = np.concatenate(
                    (self.solver_array, np.array(sol.y).T)
                )
                self.total_time = np.concatenate((self.total_time, np.array(sol.t)))

        utils.plots(self.total_time, self.solver_array.T, self.plots)

    def save_json(self):
        glide_vars = {
            "alpha_d": self.alpha_d,
            "glide_dir": self.glider_direction,
            "glide_angle_deg": self.glide_angle_deg,
            "lim1": self.lim1,
            "lim2": self.lim2,
            "theta_d": self.theta_d,
            "mb_d": self.mb_d,
            "v1_d": self.v1_d,
            "v3_d": self.v3_d,
            "m0_d": self.m0_d,
            "rp1_d": self.rp1_d,
            "rp3": self.rp3,
            "rb1": self.rb1,
            "rb3": self.rb3,
            # "rw1": self.rw1,
            # "rw2": self.rw2,
            # "rw3": self.rw3,
            "Pp1_d": self.Pp1_d,
            "Pp3_d": self.Pp3_d,
            "phi": self.phi,
            "theta0": self.theta0,
            "psi": self.psi,
            "Mf": self.Mf.tolist(),
            "M": self.M.tolist(),
            "J": self.J.tolist(),
            "KL": self.KL,
            "KL0": self.KL0,
            "KD": self.KD,
            "KD0": self.KD0,
            "KM": self.KM,
            "KM0": self.KM0,
            "KOmega1": self.KOmega1,
            "KOmega2": self.KOmega2,
            "desired_glide_speed": self.V_d,
            "ballast_rate": self.ballast_rate,
            "mh": self.mh,
            "mb": self.mb,
            "mw": self.mw,
            "mm": self.mm,
            "ms": self.ms,
            "m": self.m,
            "m0": self.m0,
            "mt": self.mt,
        }

        utils.save_json(glide_vars)

    def solve_ode(self, z0):
        def dvdt(t, y):
            eom = Dynamics(y)
            D = eom.set_eom()
            return D

        self.t = np.linspace(375 * (self.i), 375 * (self.i + 1))

        sol = solve_ivp(
            dvdt,
            t_span=(min(self.t), max(self.t)),
            y0=z0,
            method="RK45",
            t_eval=self.t,
            dense_output=False,
        )

        return sol


if __name__ == "__main__":
    Z = Vertical_Motion()
    Z.set_desired_trajectory()
