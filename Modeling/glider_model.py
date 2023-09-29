import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
import utils
import os
import json
from Parameters.slocum import SLOCUM_PARAMS
from Modeling.dynamics import Dynamics


class Vertical_Motion:
    def __init__(self):
        self.mass_params = SLOCUM_PARAMS.GLIDER_CONFIG
        self.hydro_params = SLOCUM_PARAMS.HYDRODYNAMICS
        self.vars = SLOCUM_PARAMS.VARIABLES

        self.initialization("SLOCUM")

    def initialization(self, glider_name="SLOCUM"):
        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        if glider_name == "SLOCUM":
            from Parameters.slocum import SLOCUM_PARAMS as P

        self.mass_params = P.GLIDER_CONFIG
        self.hydro_params = P.HYDRODYNAMICS
        self.vars = P.VARIABLES

        self.mh = self.mass_params.HULL_MASS
        self.mw = self.mass_params.FIXED_POINT_MASS
        self.mb = self.mass_params.BALLAST_MASS
        self.mm = self.mass_params.INT_MOVABLE_MASS

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

        [self.rp1, self.rp2, self.rp3] = [self.vars.rp1, self.vars.rp2, self.vars.rp3]
        [self.rb1, self.rb2, self.rb3] = [self.vars.rb1, self.vars.rb2, self.vars.rb3]
        [self.rw1, self.rw2, self.rw3] = [self.vars.rw1, self.vars.rw2, self.vars.rw3]

        self.glide_angle_deg = self.vars.GLIDE_ANGLE
        self.V_d = self.vars.SPEED
        self.ballast_rate = self.vars.BALLAST_RATE

        self.set_first_run_params()

    def set_first_run_params(self):
        self.n1_0 = np.array([[0.0, 0.0, 0.0]])
        self.Omega0 = np.array([[0.0, 0.0, 0.0]])
        self.phi = self.mass_params.PHI
        self.theta0 = self.mass_params.THETA
        self.psi = self.mass_params.PSI

        self.Pp = [0.0, 0.0, 0.0]
        self.Pb = [0.0, 0.0, 0.0]
        self.Pw = [0.0, 0.0, 0.0]

    def set_desired_trajectory(self):
        # sourcery skip: use-fstring-for-formatting
        self.E_i_d = np.array(
            [
                math.radians(-self.glide_angle_deg),
                math.radians(self.glide_angle_deg),
                math.radians(-self.glide_angle_deg),
                math.radians(self.glide_angle_deg),
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

        for i in range(1):
            e_i_d = self.E_i_d[i]

            print(
                "Iteration {} | Desired glide angle in deg = {}".format(
                    i, math.degrees(e_i_d)
                )
            )

            if e_i_d > 0:
                self.glider_direction = "U"
                self.ballast_rate = -abs(self.ballast_rate)
                print("Glider moving in upward direction\n")

            elif e_i_d < 0:
                self.glider_direction = "D"
                self.ballast_rate = abs(self.ballast_rate)
                print("Glider moving in downward direction\n")

            self.alpha_d = (
                (1 / 2)
                * (self.KL / self.KD)
                * math.tan(e_i_d)
                * (
                    -1
                    + math.sqrt(
                        1
                        - 4
                        * (self.KD / math.pow(self.KL, 2))
                        * (1 / math.tan(e_i_d))
                        * (self.KD0 * (1 / math.tan(e_i_d)) + self.KL0)
                    )
                )
            )

            self.mb_d = (self.m - self.mh - self.mm) + (1 / self.g) * (
                -math.sin(e_i_d) * (self.KD0 + self.KD * math.pow(self.alpha_d, 2))
                + math.cos(e_i_d) * (self.KL0 + self.KL * self.alpha_d)
            ) * math.pow(self.V_d, 2)

            self.theta_d = e_i_d + self.alpha_d

            self.v1_d = self.V_d * math.cos(self.alpha_d)
            self.v3_d = self.V_d * math.sin(self.alpha_d)

            self.Pp1_d = self.mm * self.v1_d
            self.Pp3_d = self.mm * self.v3_d

            self.Pb1_d = self.mb * self.v1_d
            self.Pb3_d = self.mb * self.v3_d

            self.m0_d = self.mb_d + self.mh + self.mm - self.m

            self.rp1_d = -self.rp3 * math.tan(self.theta_d) + (
                1 / (self.mm * self.g * math.cos(self.theta_d))
            ) * (
                (self.Mf[2, 2] - self.Mf[0, 0]) * self.v1_d * self.v3_d
                + (self.KM0 + self.KM * self.alpha_d) * math.pow(self.V_d, 2)
            )

            self.save_json()

            # These are the initial conditions at every peak of the sawtooth trajectory
            if i == 0:
                self.z_in = [
                    self.n1_0.transpose(),
                    self.Omega0.transpose(),
                    np.array([[self.v1_d, 0.0, self.v3_d]]).transpose(),
                    np.array([[self.rp1_d, self.rp2, self.rp3]]).transpose(),
                    np.array([[self.rb1, self.rb2, self.rb3]]).transpose(),
                    np.array([[self.Pp1_d, self.Pp[1], self.Pp3_d]]).transpose(),
                    np.array([[self.Pb1_d, self.Pb[1], self.Pb3_d]]).transpose(),
                    np.array([self.Pw]),
                    self.mb_d,
                    self.theta0,
                ]

            else:
                self.z_in = None
                # empty arrays

            # print(z)

            eom = Dynamics(self.z_in)
            self.Z = eom.set_eom()
            
            self.solve_ode()

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
            "Pp1_d": self.Pp1_d,
            "Pp3_d": self.Pp3_d,
            "Pb1_d": self.Pb1_d,
            "Pb3_d": self.Pb3_d,
            "m0_d": self.m0_d,
            "rp1_d": self.rp1_d,
            "rp2": self.rp2,
            "rp3": self.rp3,
            "rb1": self.rb1,
            "rb2": self.rb2,
            "rb3": self.rb3,
            "rw1": self.rw1,
            "rw2": self.rw2,
            "rw3": self.rw3,
            "phi": self.phi,
            "theta0": self.theta0,
            "psi": self.psi,
            "Mf": self.Mf.tolist(),
        }

        glide_vars["M"] = self.M.tolist()
        glide_vars["J"] = self.J.tolist()
        glide_vars["KL"] = self.KL
        glide_vars["KL0"] = self.KL0
        glide_vars["KD"] = self.KD
        glide_vars["KD0"] = self.KD0
        glide_vars["KM"] = self.KM
        glide_vars["KM0"] = self.KM0
        glide_vars["KOmega1"] = self.KOmega1
        glide_vars["KOmega2"] = self.KOmega2
        glide_vars["desired_glide_speed"] = self.V_d
        glide_vars["ballast_rate"] = self.ballast_rate
        glide_vars["mh"] = self.mh
        glide_vars["mb"] = self.mb
        glide_vars["mw"] = self.mw
        glide_vars["mm"] = self.mm
        glide_vars["m"] = self.m
        glide_vars["m0"] = self.m0
        glide_vars["mt"] = self.mt

        utils.save_json(glide_vars)

    def solve_ode(self):
        import pdb

        pdb.set_trace()
        t = np.linspace(0, 200, 500)
        sol = solve_ivp(self.Z, t_span=(0, max(t)), y0=self.z_in, t_eval=t)
        print(sol)


if __name__ == "__main__":
    Z = Vertical_Motion()
    Z.set_desired_trajectory()
    # Z.save_json()
