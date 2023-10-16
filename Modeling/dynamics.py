import numpy as np
import math
import utils
from Parameters.slocum import SLOCUM_PARAMS


class Dynamics:
    def __init__(self, z):
        self.initialization()

        self.n1 = np.array([z[0:3]]).T
        self.Omega = np.array([z[3:6]]).T
        self.v = np.array([z[6:9]]).T
        self.rp = np.array([z[9:12]]).T
        self.rb = np.array([z[12:15]]).T
        self.Pp = np.array([z[15:18]]).T
        self.Pb = np.array([z[18:21]]).T
        self.mb = z[21]
        self.theta = z[24]
        print(math.sqrt(math.pow(self.v[0][0], 2) + math.pow(self.v[2][0], 2)))

        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        self.controls = SLOCUM_PARAMS.CONTROLS

        # if self.glide_dir == "U":
        #     if self.mb <= self.mb_d:
        #         self.mb = self.mb_d

        # elif self.glide_dir == "D":
        #     if self.mb >= self.mb_d:
        #         self.mb = self.mb_d

        # if self.glide_dir == "U":
        #     if self.rp[0] <= self.rp1_d:
        #         self.rp[0] = self.rp1_d

        # elif self.glide_dir == "D":
        #     if self.rp[0] >= self.rp1_d:
        #         self.rp[0] = self.rp1_d

        self.rp_c = utils.unit_vecs(self.rp)
        self.rb_c = utils.unit_vecs(self.rb)

        

        self.R, self.R_T = self.transformation()

        self.set_force_torque()
        
        self.set_limits()

        # if self.glide_dir == "U":
        #     self.rp_dot = (
        #         (1 / self.mm) * self.Pp - self.v - np.cross(self.Omega, self.rp, axis=0)
        #     )

        #     if self.rp[0] <= self.rp1_d:
        #         self.rp_dot = np.array([self.Z3]).transpose()

        # elif self.glide_dir == "D":
        #     self.rp_dot = (
        #         (1 / self.mm) * self.Pp + self.v - np.cross(self.Omega, self.rp, axis=0)
        #     )

        #     if self.rp[0] >= self.rp1_d:
        #         self.rp_dot = np.array([self.Z3]).transpose()


        # self.rb_dot = np.array([self.Z3]).transpose()

        # self.Pb = self.mb * (
        #     self.v + np.cross(self.Omega, self.rb, axis=0) + self.rb_dot
        # )

        self.m0 = self.mh + self.mw + self.mb + self.mm - self.m

        # if self.glide_dir == "U":
        #     if self.mb <= self.mb_d:
        #         self.ballast_rate = 0
        # elif self.glide_dir == "D":
        #     if self.mb >= self.mb_d:
        #         self.ballast_rate = 0

    def initialization(self):
        var = utils.load_json("vars/glider_variables.json")

        self.alpha_d = var["alpha_d"]
        self.glide_dir = var["glide_dir"]
        self.glide_angle_deg = var["glide_angle_deg"]
        self.lim1 = var["lim1"]
        self.lim2 = var["lim2"]
        self.theta_d = var["theta_d"]
        self.mb_d = var["mb_d"]
        self.v1_d = var["v1_d"]
        self.v3_d = var["v3_d"]
        self.m0_d = var["m0_d"]
        self.rp1_d = var["rp1_d"]
        self.rp3 = var["rp3"]
        self.rb1 = var["rb1"]
        self.rb3 = var["rb3"]
        # self.rw1 = var["rw1"]
        # self.rw2 = var["rw2"]
        # self.rw3 = var["rw3"]
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
        self.ms = var["ms"]
        self.m = var["m"]
        self.m0 = var["m0"]
        self.mt = var["mt"]

    def transformation(self):
        # n2 = np.array([[self.phi, self.theta, self.psi]]).transpose()
        # eta = np.array(
        #     [self.n1.transpose(), n2.transpose()], dtype=np.float32
        # ).transpose()
        # self.v = np.array([self.v, self.Omega], dtype=np.float32)
        # tau1 = np.array([[X, Y, Z]]).transpose()
        # tau2 = np.array([[K, M, N]]).transpose()
        # tau = np.array([tau1.transpose(), tau2.transpose()],
        #                dtype=np.float32).transpose()

        J1_n2, J2_n2 = utils.transformationMatrix(self.phi, self.theta, self.psi)

        # n1_dot = J1_n2 * self.v
        # n2_dot = J2_n2 * self.Omega

        # self.eta_dot = np.array([n1_dot, n2_dot]).transpose()

        # self.eta_dot_T = np.linalg.inv(self.R)

        # Replace R with J1_n2 later

        R = J1_n2
        # self.R_T = np.linalg.inv(self.R)
        R_T = R.T

        return R, R_T

    def set_force_torque(self):
        self.alpha = math.atan(self.v[2][0] / self.v[0][0])

        L = (self.KL0 + self.KL * self.alpha) * (
            math.pow(self.v[0][0], 2) + math.pow(self.v[2][0], 2)
        )
        D = (self.KD0 + self.KD * (math.pow(self.alpha, 2))) * (
            math.pow(self.v[0][0], 2) + math.pow(self.v[2][0], 2)
        )
        MDL = (self.KM0 + self.KM * self.alpha) * (
            math.pow(self.v[0][0], 2) + math.pow(self.v[2][0], 2)
        )
        +self.KOmega1 * self.Omega[1][0]
        +self.KOmega2 * math.pow(self.Omega[1][0], 2)

        self.F_ext = np.array([[-D, 0, -L]]).transpose()  # same as X, Y, Z

        self.T_ext = np.array([[0, MDL, 0]]).transpose()  # same as K, M, N

    def set_limits(self):
        if self.glide_dir == "D":
            self.rp_dot = (
                (1 / self.mm) * self.Pp
                + self.v
                - np.cross(self.Omega, self.rp, axis=0)
            )

            if self.rp[0] >= self.rp1_d:
                self.rp[0] = self.rp1_d
                self.rp_dot = np.array([self.Z3]).transpose()

            if self.mb >= self.mb_d:
                self.mb = self.mb_d
                self.ballast_rate = 0

            # Write for mw also

            # self.rb_dot = (1/mb) * self.Pb + self.v - np.cross(self.Omega, self.rb) # Uncomment if ballast mass moves
            # Assume ballast mass does not move
            self.rb_dot = np.array([self.Z3]).transpose()

            self.rw_dot = np.array([self.Z3]).transpose()

        elif self.glide_dir == "U":
            self.rp_dot = (
                (1 / self.mm) * self.Pp
                - self.v
                - np.cross(self.Omega, self.rp, axis=0)
            )

            if self.rp[0] <= self.rp1_d:
                self.rp[0] = self.rp1_d
                self.rp_dot = np.array([self.Z3]).transpose()

            if self.mb <= self.mb_d:
                self.mb = self.mb_d
                self.ballast_rate = 0

            # Write for mw also

            # self.rb_dot = (1/mb) * self.Pb - self.v - np.cross(self.Omega, self.rb) # Uncomment if ballast mass moves
            # Assume ballast mass does not move
            self.rb_dot = np.array([self.Z3]).transpose()

            self.rw_dot = np.array([self.Z3]).transpose()

        self.Pb = self.mb * (
            self.v + np.cross(self.Omega, self.rb, axis=0) + self.rb_dot
        )

    def control_transformation(self):
        if self.glide_dir == "D":
            self.wp = np.array(
                [[self.controls.wp1, self.controls.wp2, self.controls.wp3]]
            ).transpose()
        elif self.glide_dir == "U":
            self.wp = np.array(
                [[-self.controls.wp1, -self.controls.wp2, -self.controls.wp3]]
            ).transpose()

        self.wb = np.array([[0, 0, 0]]).transpose()

        self.ww = np.array([[0, 0, 0]]).transpose()

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
                np.cross(self.M @ (self.v) + self.Pp + self.Pb, self.Omega, axis=0)
                + self.m0 * self.g * np.matmul(self.R_T, self.k_hat)
                + self.F_ext
            )
            - np.cross(self.Omega, self.rp_dot, axis=0)
            - np.cross(
                np.linalg.inv(self.J)
                @ (
                    np.cross(
                        np.matmul(self.J, self.Omega)
                        + np.matmul(self.rp_c, self.Pp)
                        + np.matmul(self.rb_c, self.Pb),
                        self.Omega,
                        axis=0,
                    )
                    + np.cross(np.matmul(self.M, self.v), self.v, axis=0)
                    + self.T_ext
                    + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
                    + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
                    + (self.mm * self.rp_c + self.mb * self.rb_c)
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
                np.cross(self.M @ (self.v) + self.Pp + self.Pb, self.Omega, axis=0)
                + self.m0 * self.g * np.matmul(self.R_T, self.k_hat)
                + self.F_ext
            )
            - np.cross(self.Omega, self.rp_dot, axis=0)
            - np.cross(
                np.linalg.inv(self.J)
                @ (
                    np.cross(
                        np.matmul(self.J, self.Omega)
                        + np.matmul(self.rp_c, self.Pp)
                        + np.matmul(self.rb_c, self.Pb),
                        self.Omega,
                        axis=0,
                    )
                    + np.cross(np.matmul(self.M, self.v), self.v, axis=0)
                    + self.T_ext
                    + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
                    + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
                    + (self.mm * self.rp_c + self.mb * self.rb_c)
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
        #         np.cross(self.M @ (self.v) + self.Pp + self.Pb + self.Pw, self.Omega, axis=0)
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
        #         + np.cross(np.matmul(self.M, self.v), self.v, axis=0)
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
        self.u = self.H @ np.array(
            [(-Zp + self.wp).flatten(), (-Zb + self.wb).flatten()]
        ).reshape(6, 1)
        self.u = self.u.reshape(2, 3)

        self.u_bar = np.array([self.u[0]]).T

        self.u_b = np.array([self.Z3]).T

        # self.u_w = np.array([self.Z3]).T

    def set_controls(self):
        self.control_transformation()

        if self.glide_dir == "D":
            if self.rp[0] >= self.rp1_d:
                self.u_bar = np.array([self.Z3]).T
            # else:
            #     self.u_bar = np.array([[0.18, 0, 0]]).T

            # if self.rs[0] >= self.rs1_d:
            #     self.u_s = np.array([self.Z3]).T
            # else:
            #     self.u_s = np.array([[0.18, 0, 0]]).T

            # if self.ms == 0:
            #     self.u_s = np.array([[0, 0, 0]]).T

        if self.glide_dir == "U":
            if self.rp[0] <= self.rp1_d:
                self.u_bar = np.array([self.Z3]).T
            # else:
            #     self.u_bar = np.array([[-0.18, 0, 0]]).T

            # if self.rs[0] <= self.rs1_d:
            #     self.u_s = np.array([self.Z3]).T
            # else:
            #     self.u_s = np.array([[-0.18, 0, 0]]).T

            # if self.ms == 0:
            #     self.u_s = np.array([[0, 0, 0]]).T

        # Write for u_b and u_w as well

        self.u_b = np.array([self.Z3]).T

    def set_eom(self):
        self.set_controls()

        T_bar = (
            np.cross(
                (self.J @ self.Omega + self.rp_c @ self.Pp + self.rb_c @ self.Pb),
                self.Omega,
                axis=0,
            )
            + np.cross(self.M @ (self.v), self.v, axis=0)
            + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
            + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
            + (self.mm * self.rp_c + self.mb * self.rb_c)
            * self.g
            @ np.matmul(self.R_T, self.k_hat)
            + self.T_ext
            - self.rp_c @ (self.u_bar)
            - (self.rb_c @ (self.u_b))
        )

        F_bar = (
            np.cross((self.M @ self.v + self.Pp + self.Pb), self.Omega, axis=0)
            + self.m0 * self.g * self.R_T @ self.k_hat
            + self.F_ext
            - self.u_bar
            - (self.u_b)
        )

        # R_dot

        n1_dot = self.R @ self.v

        Omega_dot = np.linalg.inv(self.J) @ T_bar

        v1_dot = np.linalg.inv(self.M) @ F_bar

        # rp_dot solved earlier

        # rb_dot solved earlier

        # rw_dot solved earlier

        # rs_dot solved earlier

        Pp_dot = self.u_bar

        Pb_dot = self.u_b

        # Ps_dot = self.u_s

        # Pw_dot = self.u_w

        mb_dot = np.array([self.ballast_rate])

        theta_dot = np.array([self.Omega[1][0]])

        # breakpoint()

        return np.array(
            [
                n1_dot.flatten(),
                Omega_dot.flatten(),
                v1_dot.flatten(),
                self.rp_dot.flatten(),
                #  self.rs_dot.flatten(),
                self.rb_dot.flatten(),
                Pp_dot.flatten(),
                #  Ps_dot.flatten(),
                Pb_dot.flatten(),
                [mb_dot[0], 0, 0],
                [theta_dot[0], 0, 0],
            ]
        ).flatten()


if __name__ == "__main__":
    equations = Dynamics()
    # equations.set_eom()
    # equations.initialization()
