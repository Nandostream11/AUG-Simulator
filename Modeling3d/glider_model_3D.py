import numpy as np
import math
from scipy.integrate import solve_ivp
import utils
from Modeling3d.dynamics_3D import Dynamics


class ThreeD_Motion:
    def __init__(self, args):
        self.wp = []
        self.args = args
        self.mode = self.args.mode
        if self.mode == "3D":
            self.cycles = 1
            self.rudder = self.args.rudder
            self.rudder_angle = math.radians(self.args.setrudder)
        else:
            self.cycles = self.args.cycle
        self.glider_name = self.args.glider
        self.info = self.args.info
        self.pid_control = self.args.pid
        self.plots = self.args.plot

        self.initialization()

        self.solver_array = np.array([])
        self.total_time = np.array([])

    def initialization(self):
        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        if self.glider_name == ("slocum") and self.mode == "3D":
            from Parameters.slocum3D import SLOCUM_PARAMS as P
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
        self.K_beta = self.hydro_params.K_beta
        self.KM = self.hydro_params.KM
        self.KM0 = self.hydro_params.KM0
        self.K_MY = self.hydro_params.K_MY
        self.K_MR = self.hydro_params.K_MR
        self.KOmega11 = self.hydro_params.KOmega11
        self.KOmega12 = self.hydro_params.KOmega12
        self.KOmega13 = self.hydro_params.KOmega13
        self.KOmega21 = self.hydro_params.KOmega21
        self.KOmega22 = self.hydro_params.KOmega22
        self.KOmega23 = self.hydro_params.KOmega23

        if self.rudder == "enable":
            self.rp2_d = 0.0
        if self.rudder == "disable":
            self.rp2_d = self.vars.rp2

        self.rp3 = self.vars.rp3
        self.rb1 = self.vars.rb1
        self.rb2 = self.vars.rb2
        self.rb3 = self.vars.rb3

        self.glide_angle_deg = self.vars.GLIDE_ANGLE
        self.V_d = self.vars.SPEED
        self.ballast_rate = self.vars.BALLAST_RATE

        self.set_first_run_params()

    def set_first_run_params(self):
        self.phi0 = math.radians(self.vars.PHI)
        self.theta0 = -math.radians(self.vars.THETA)
        self.psi0 = math.radians(self.vars.PSI)
        self.Omega0 = [0.0046, 0.0025, 0.0077]

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

            self.beta_d = math.radians(self.vars.BETA)

            self.mb_d = (self.m - self.mh - self.mm) + (1 / self.g) * (
                -math.sin(self.e_i_d) * (self.KD0 + self.KD * math.pow(self.alpha_d, 2))
                + math.cos(self.e_i_d) * (self.KL0 + self.KL * self.alpha_d)
            ) * math.pow(self.V_d, 2)

            self.m0_d = self.mb_d + self.mh + self.mm - self.m

            self.theta_d = self.e_i_d + self.alpha_d

            self.v1_d = self.V_d * math.cos(self.alpha_d) * math.cos(self.beta_d)
            self.v2_d = self.V_d * math.sin(self.beta_d)
            self.v3_d = self.V_d * math.sin(self.alpha_d) * math.cos(self.beta_d)

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
                    "Desired longitudinal position of internal movable mass in cm = {}".format(
                        self.rp1_d * 100
                    )
                )

            self.save_json()

            # Initial conditions for spiral motion

            if i == 0:
                self.z_in = np.concatenate(
                    [
                        [0.0, 0.0, 0.0],
                        [self.Omega0[0], self.Omega0[1], self.Omega0[2]],
                        [self.v1_d, self.v2_d, self.v3_d],
                        [0.0, 0.0, self.rp3],
                        [self.rb1, self.rb2, self.rb3],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [self.mb_d, 0, 0],
                        [self.phi0, self.theta0, self.psi0],
                    ]
                ).ravel()

            else:
                self.z_in = self.solver_array[-1]

            self.t = np.linspace(2000 * (i), 2000 * (i + 1), 1000)

            sol, w = self.solve_ode(self.z_in, self.t)

            if i == 0:
                self.solver_array = sol.y.T
                self.total_time = sol.t
                self.wp = w
            else:
                self.solver_array = np.concatenate((self.solver_array, sol.y.T))
                self.total_time = np.concatenate((self.total_time, sol.t))
                self.wp = np.concatenate((self.wp, w))

            if self.mode == "3D":
                v = math.sqrt(
                    math.pow(self.solver_array[-1][6], 2)
                    + math.pow(self.solver_array[-1][7], 2)
                    + math.pow(self.solver_array[-1][8], 2)
                )
                beta = math.asin(self.solver_array[-1][7] / v)
                alpha = math.atan(self.solver_array[-1][8] / self.solver_array[-1][6])
                R = (
                    v
                    * math.cos(self.solver_array[-1][-2] - alpha)
                    / self.solver_array[-1][5]
                )
                print(
                    "\nEquilibrium roll angle of glider: {} deg".format(
                        math.degrees(self.solver_array[-1][-3])
                    )
                )
                print(
                    "Equilibrium pitch angle of glider: {} deg".format(
                        math.degrees(self.solver_array[-1][-2])
                    )
                )
                print("Sideslip angle of glider: {} deg".format(math.degrees(beta)))
                print("Equilibrium glide speed: {} m/s".format(v))
                print("Radius : {} m".format(R))

        utils.plots(self.total_time, self.solver_array.T, self.plots)

    def save_json(self):
        glide_vars = {
            "glide_dir": self.glider_direction,
            "glide_angle_deg": self.glide_angle_deg,
            "lim1": self.lim1,
            "lim2": self.lim2,
            "alpha_d": self.alpha_d,
            "beta_d": self.beta_d,
            "theta_d": self.theta_d,
            "mb_d": self.mb_d,
            "v1_d": self.v1_d,
            "v2_d": self.v2_d,
            "v3_d": self.v3_d,
            "m0_d": self.m0_d,
            "rp1_d": self.rp1_d,
            "rp2": self.rp2_d,
            "rp3": self.rp3,
            "rb1": self.rb1,
            "rb2": self.rb2,
            "rb3": self.rb3,
            "phi0": self.phi0,
            "theta0": self.theta0,
            "psi0": self.psi0,
            "Mf": self.Mf.tolist(),
            "M": self.M.tolist(),
            "J": self.J.tolist(),
            "KL": self.KL,
            "KL0": self.KL0,
            "KD": self.KD,
            "KD0": self.KD0,
            "K_beta": self.K_beta,
            "KM": self.KM,
            "KM0": self.KM0,
            "K_MY": self.K_MY,
            "K_MR": self.K_MR,
            "KOmega11": self.KOmega11,
            "KOmega12": self.KOmega12,
            "KOmega13": self.KOmega13,
            "KOmega21": self.KOmega21,
            "KOmega22": self.KOmega22,
            "KOmega23": self.KOmega23,
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
            "pid_control": self.pid_control,
            "rudder": self.rudder,
            "rudder_angle": self.rudder_angle,
        }

        utils.save_json(glide_vars, "vars/3d_glider_variables.json")

    def solve_ode(self, z0, time):
        def dvdt(t, y):
            global inner_func

            def inner_func(t, y):
                eom = Dynamics(y)
                D = eom.set_eom()
                return D

            Dr = inner_func(t, y)
            return Dr[:-3]

        sol = solve_ivp(
            dvdt,
            t_span=(min(time), max(time)),
            y0=z0,
            method="RK45",
            t_eval=time,
            dense_output=False,
        )

        w = np.array([inner_func(time[i], sol.y.T[i, :]) for i in range(len(time))])[
            :, -2
        ]

        return sol, w


if __name__ == "__main__":
    Z = ThreeD_Motion()
    Z.set_desired_trajectory()
