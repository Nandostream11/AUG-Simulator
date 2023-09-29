import numpy as np
import math
import utils
from Parameters.slocum import SLOCUM_PARAMS
import json

# f_gravity = mt * g * np.dot(J1_n2.transpose(), k_hat)
# f_buoyancy = -m * g * np.dot(J1_n2.transpose(), k_hat)
# tau_gravity = np.cross(r_cg, mt * g * np.dot(J1_n2.transpose, k_hat))
# tau_buoyancy = 0

# Terms to be defined:
# mb_d
# rp_d
# ballast_rate
# phi, theta, psi -- all angles in degrees
# masses


class Dynamics:
    def __init__(self, z):
        # b same as n1, v same as v1, ohmega same as v2

        self.n1 = z[0]
        self.Omega = z[1]
        self.v1 = z[2]
        self.rp = z[3]
        self.rb = z[4]
        self.rw = z[5]
        self.Pp = z[6]
        self.Pb = z[7]
        self.Pw = z[8]
        self.mb = z[9][0]
        self.theta = z[10][0]

        self.rp_c = utils.unit_vecs(self.rp)
        self.rb_c = utils.unit_vecs(self.rb)
        self.rw_c = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])

        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        self.controls = SLOCUM_PARAMS.CONTROLS

        self.initialization()

        self.set_limits()

        self.transformation()

        self.set_force_torque()

    def initialization(self):
        var = utils.load_json("vars/glider_variables.json")

        self.alpha_d = var["alpha_d"]
        self.glide_dir = var["glide_dir"]
        self.glide_angle_deg = var["glide_dir"]
        self.lim1 = var["lim1"]
        self.lim2 = var["lim2"]
        self.theta_d = var["theta_d"]
        self.mb_d = var["mb_d"]
        self.v1_d = var["v1_d"]
        self.v3_d = var["v3_d"]
        self.Pp1_d = var["Pp1_d"]
        self.Pp3_d = var["Pp3_d"]
        self.m0_d = var["m0_d"]
        self.rp1_d = var["rp1_d"]
        self.rp2 = var["rp2"]
        self.rp3 = var["rp3"]
        self.rb1 = var["rb1"]
        self.rb2 = var["rb2"]
        self.rb3 = var["rb3"]
        self.rw1 = var["rw1"]
        self.rw2 = var["rw2"]
        self.rw3 = var["rw3"]
        self.phi = var["phi"]
        self.theta0 = var["theta0"]
        self.psi = var["psi"]

        self.Mf = np.array(var["Mf"])
        self.M = np.array(var["M"])
        self.J = np.array(var["J"])
        self.KL = var["KL"]
        self.KL0 = var["KL0"]
        self.KM = var["KM"]
        self.KM0 = var["KM0"]
        self.KD = var["KD"]
        self.KD0 = var["KD0"]
        self.KOmega1 = var["KOmega1"]
        self.KOmega2 = var["KOmega2"]
        self.V_d = var["desired_glide_speed"]
        self.ballast_rate = var["ballast_rate"]
        self.mh = var["mh"]
        self.mb = var["mb"]
        self.mw = var["mw"]
        self.mm = var["mm"]
        self.m = var["m"]
        self.m0 = var["m0"]
        self.mt = var["mt"]

    def transformation(self):
        n2 = np.array([[self.phi, self.theta, self.psi]]).transpose()
        eta = np.array(
            [self.n1.transpose(), n2.transpose()], dtype=np.float32
        ).transpose()
        self.v = np.array([self.v1, self.Omega], dtype=np.float32)
        # tau1 = np.array([[X, Y, Z]]).transpose()
        # tau2 = np.array([[K, M, N]]).transpose()
        # tau = np.array([tau1.transpose(), tau2.transpose()],
        #                dtype=np.float32).transpose()

        J1_n2, J2_n2 = utils.transformationMatrix(self.phi, self.theta, self.psi)

        # n1_dot = J1_n2 * self.v1
        # n2_dot = J2_n2 * self.Omega

        # self.eta_dot = np.array([n1_dot, n2_dot]).transpose()

        # self.eta_dot_T = np.linalg.inv(self.R)

        # Replace R with J1_n2 later

        self.R = J1_n2
        self.R_T = np.linalg.inv(self.R)

    def set_force_torque(self):
        self.alpha = np.arctan(self.v1[2][0] / self.v1[0][0])

        L = (self.KL0 + self.KL * self.alpha) * (
            math.pow(self.v1[0][0], 2) + math.pow(self.v1[2][0], 2)
        )
        D = (self.KD0 + self.KD * (math.pow(self.alpha, 2))) * (
            math.pow(self.v1[0][0], 2) + math.pow(self.v1[2][0], 2)
        )
        MDL = (self.KM0 + self.KM * self.alpha) * (
            math.pow(self.v1[0][0], 2)
            + math.pow(self.v1[2][0], 2)
            + self.KOmega1 * self.Omega[1][0]
            + self.KOmega2 * math.pow(self.Omega[1][0], 2)
        )

        self.F_ext = np.array([[-D, 0, -L]]).transpose()  # same as X, Y, Z

        self.T_ext = np.array([[0, MDL, 0]]).transpose()  # same as K, M, N

        return self.F_ext, self.T_ext

    def set_limits(self):
        if self.glide_dir == "D":
            if self.rp[0] >= self.rp1_d:
                self.rp[0] = self.rp1_d

            if self.mb >= self.mb_d:
                self.mb = self.mb_d
                self.ballast_rate = 0

            # Write for mw also

            self.rp_dot = (
                (1 / self.mm) * self.Pp
                + self.v1
                - np.cross(self.Omega, self.rp, axis=0)
            )

            # self.rb_dot = (1/mb) * self.Pb + self.v1 - np.cross(self.Omega, self.rb) # Uncomment if ballast mass moves
            # Assume ballast mass does not move
            self.rb_dot = np.array([self.Z3]).transpose()

            self.rw_dot = np.array([self.Z3]).transpose()

        elif self.glide_dir == "U":
            if self.rp[0] <= self.rp1_d:
                self.rp[0] = self.rp1_d

            if self.mb <= self.mb_d:
                self.mb = self.mb_d
                self.ballast_rate = 0

            # Write for mw also

            self.rp_dot = (
                (1 / self.mm) * self.Pp
                - self.v1
                - np.cross(self.Omega, self.rp, axis=0)
            )

            # self.rb_dot = (1/mb) * self.Pb - self.v1 - np.cross(self.Omega, self.rb) # Uncomment if ballast mass moves
            # Assume ballast mass does not move
            self.rb_dot = np.array([self.Z3]).transpose()

            self.rw_dot = np.array([self.Z3]).transpose()

    def control_transformation(self):
        if self.glide_dir == "D":
            self.wp = np.array(
                [[self.controls.wp1, self.controls.wp2, self.controls.wp3]]
            ).transpose()
        elif self.glide_dir == "U":
            self.wp = np.array(
                [[-self.controls.wp1, -self.controls.wp2, -self.controls.wp3]]
            ).transpose()

        self.wb = np.array(
            [[self.controls.wb1, self.controls.wb2, self.controls.wb3]]
        ).transpose()

        self.ww = np.array(
            [[self.controls.ww1, self.controls.ww2, self.controls.ww3]]
        ).transpose()

        # 3x3 matrix when mw not equal to 0 -- needs to be updated
        # self.F = np.array(
        #     [
        #         [
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rp_c
        #             + (1 / self.mm) * self.I3,
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rb_c,
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rw_c,
        #         ],
        #         [
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rp_c,
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rb_c
        #             + (1 / self.mb) * self.I3,
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rw_c,
        #         ],
        #         [
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rp_c,
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rb_c,
        #             np.linalg.inv(self.M)
        #             - self.rp_c * np.linalg.inv(self.J) * self.rw_c
        #             + (1 / self.mw) * self.I3,
        #         ],
        #     ]
        # )

        # 2x2 matrix when mw = 0
        self.F = np.array(
            [
                [
                    np.linalg.inv(self.M)
                    - self.rp_c @ (np.linalg.inv(self.J)) @ (self.rp_c)
                    + (1 / self.mm) * self.I3,
                    np.linalg.inv(self.M)
                    - self.rp_c @ (np.linalg.inv(self.J)) @ (self.rb_c),
                ],
                [
                    np.linalg.inv(self.M)
                    - self.rb_c @ (np.linalg.inv(self.J)) @ (self.rp_c),
                    np.linalg.inv(self.M)
                    - self.rb_c @ (np.linalg.inv(self.J)) @ (self.rb_c)
                    + (1 / self.mb) * self.I3,
                ],
            ]
        )

        # self.F = np.linalg.det(self.F)  # But why?

        self.H = np.linalg.inv(self.F)

        Zp = (
            -np.linalg.inv(self.M)
            @ (
                np.cross(
                    self.M @ (self.v1) + self.Pp + self.Pb + self.Pw, self.Omega, axis=0
                )
                + self.m0 * self.g * np.matmul(self.R_T, self.k_hat)
                + self.F_ext
            )
            - np.cross(self.Omega, self.rp_dot, axis=0)
            - np.linalg.inv(self.J)
            @ (
                np.cross(
                    np.matmul(self.J, self.Omega)
                    + np.matmul(self.rp_c, self.Pp)
                    + np.matmul(self.rb_c, self.Pb)
                    + np.matmul(self.rw_c, self.Pw),
                    self.Omega,
                    axis=0,
                )
                + np.cross(np.matmul(self.M, self.v1), self.v1, axis=0)
                + self.T_ext
            )
            + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
            + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
            + np.cross(np.cross(self.Omega, self.rw, axis=0), self.Pw, axis=0)
            + np.cross(
                (
                    (self.mm * self.rp_c + self.mb * self.rb_c + self.mw * self.rw_c)
                    * self.g
                    @ (self.R_T @ self.k_hat)
                ),
                self.rp,
                axis=0,
            )
        )

        Zb = (
            -np.linalg.inv(self.M)
            @ (
                np.cross(
                    self.M @ (self.v1) + self.Pp + self.Pb + self.Pw, self.Omega, axis=0
                )
                + self.m0 * self.g * np.matmul(self.R_T, self.k_hat)
                + self.F_ext
            )
            - np.cross(self.Omega, self.rp_dot, axis=0)
            - np.linalg.inv(self.J)
            @ (
                np.cross(
                    np.matmul(self.J, self.Omega)
                    + np.matmul(self.rp_c, self.Pp)
                    + np.matmul(self.rb_c, self.Pb)
                    + np.matmul(self.rw_c, self.Pw),
                    self.Omega,
                    axis=0,
                )
                + np.cross(np.matmul(self.M, self.v1), self.v1, axis=0)
                + self.T_ext
            )
            + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
            + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
            + np.cross(np.cross(self.Omega, self.rw, axis=0), self.Pw, axis=0)
            + np.cross(
                (
                    (self.mm * self.rp_c + self.mb * self.rb_c + self.mw * self.rw_c)
                    * self.g
                    @ (self.R_T @ self.k_hat)
                ),
                self.rb,
                axis=0,
            )
        )

        # Uncomment below if mw is not equal to 0
        # Zw = (
        #     -np.linalg.inv(self.M)
        #     @ (
        #         np.cross(self.M @ (self.v1) + self.Pp + self.Pb + self.Pw, self.Omega, axis=0)
        #         + self.m0 * self.g * np.matmul(self.R_T, self.k_hat)
        #         + self.F_ext
        #     )
        #     - np.cross(self.Omega, self.rp_dot, axis=0)
        #     - np.linalg.inv(self.J)
        #     @ (
        #         np.cross(
        #             np.matmul(self.J, self.Omega)
        #             + np.matmul(self.rp_c, self.Pp)
        #             + np.matmul(self.rb_c, self.Pb)
        #             + np.matmul(self.rw_c, self.Pw),
        #             self.Omega, axis=0
        #         )
        #         + np.cross(np.matmul(self.M, self.v1), self.v1, axis=0)
        #         + self.T_ext
        #     )
        #     + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
        #     + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
        #     + np.cross(np.cross(self.Omega, self.rw, axis=0), self.Pw, axis=0)
        #     + np.cross(
        #         (
        #             (self.mm * self.rp_c + self.mb * self.rb_c + self.mw * self.rw_c)
        #             * self.g
        #             @ (self.R_T @ self.k_hat)
        #         ),
        #         self.rw, axis=0
        #     )
        # )

        self.H = np.concatenate(
            (
                np.concatenate((self.H[0][0], self.H[0][1]), axis=1),
                np.concatenate((self.H[1][0], self.H[1][1]), axis=1),
            ),
            axis=0,
        )

        # self.u = self.H * np.array([[-Zp + self.wp], [-Zb + self.wb], [-Zw + self.ww]]) # uncomment if mw is not 0
        self.u = self.H @ np.array([(-Zp + self.wp), (-Zb + self.wb)]).reshape(6, 1)
        self.u = self.u.reshape(2, 3)

        self.u_bar = np.array([self.u[0]]).T

        self.u_b = np.array([self.Z3]).T
        self.u_w = np.array([self.Z3]).T

    def set_controls(self):
        self.control_transformation()

        if self.glide_dir == "D":
            if self.rp[0] >= self.rp1_d:
                self.u_bar = np.array([self.Z3]).T
            else:
                self.u_bar = self.u[0]

        if self.glide_dir == "U":
            if self.rp[0] <= self.rp1_d:
                self.u_bar = np.array([self.Z3]).T
            else:
                self.u_bar = self.u[0]

        # Write for u_b and u_w as well

    def set_eom(self):
        self.set_controls()

        T_bar = (
            np.cross(
                (
                    self.J @ self.Omega
                    + self.rp_c @ self.Pp
                    + self.rb_c @ self.Pb
                    + self.rw_c @ self.Pw
                ),
                self.Omega,
                axis=0,
            )
            + np.cross(self.M @ (self.v1), self.v1, axis=0)
            + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
            + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
            + np.cross(np.cross(self.Omega, self.rw, axis=0), self.Pw, axis=0)
            + (self.mm * self.rp_c + self.mb * self.rb_c + self.mw * self.rw_c)
            * self.g
            @ np.matmul(self.R_T, self.k_hat)
            + self.T_ext
            - self.rp_c @ (self.u_bar)
            - (self.rb_c @ (self.u_b) + self.rw_c @ (self.u_w))
        )

        F_bar = (
            np.cross(
                (self.M @ self.v1 + self.Pp + self.Pb + self.Pw), self.Omega, axis=0
            )
            + self.m0 * self.g * self.R_T @ self.k_hat
            + self.F_ext
            - self.u_bar
            - (self.u_b + self.u_w)
        )

        # R_dot

        n1_dot = self.R @ self.v1

        Omega_dot = np.linalg.inv(self.J) @ T_bar

        v1_dot = np.linalg.inv(self.M) @ F_bar

        # rp_dot solved earlier

        # rb_dot solved earlier

        # rw_dot solved earlier

        Pp_dot = self.u_bar

        Pb_dot = self.u_b

        Pw_dot = self.u_w

        mb_dot = np.array([self.ballast_rate])

        theta_dot = np.array([self.Omega[1][0]])

        return [
            n1_dot,
            Omega_dot,
            v1_dot,
            self.rp_dot,
            self.rb_dot,
            self.rw_dot,
            Pp_dot,
            Pb_dot,
            Pw_dot,
            mb_dot,
            theta_dot,
        ]


if __name__ == "__main__":
    equations = Dynamics()
    # equations.set_eom()
    equations.initialization()
