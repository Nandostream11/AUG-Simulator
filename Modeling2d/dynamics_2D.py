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
        # self.Pp = np.array([z[15:18]]).T
        # self.Pb = np.array([z[18:21]]).T
        self.rp_dot = np.array([z[15:18]]).T
        self.rb_dot = np.array([z[18:21]]).T
        self.mb = z[21]
        self.theta = z[25]
        self.w1 = 0.0  # <-- Add this line before the PID block

        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()
        
        print(math.degrees(self.theta))

        if self.pid_control == "enable":
            self.w1, self.theta_prev, self.error = utils.PID(
                0.05,
                0.0,
                0.0005,
                self.theta_d,
                self.theta,
                self.theta_prev,
                0.1,
                -self.Omega[1],
            )
            
            # breakpoint()
            # print(self.rp1)

            pid_vars = {
                "theta_prev": self.theta_prev.tolist(),
            }
            utils.save_json(pid_vars, "vars/pid_variables.json")

        else:
            self.rp1 = 0.0

        self.controls = SLOCUM_PARAMS.CONTROLS

        self.rp_c = utils.unit_vecs(self.rp)
        self.rb_c = utils.unit_vecs(self.rb)

        self.set_force_torque()

        self.R, self.R_T = self.transformation()

        self.m0 = self.mh + self.mw + self.mb + self.mm - self.m

        if self.glide_dir == "U":
            if self.mb <= self.mb_d:
                self.ballast_rate = 0
        elif self.glide_dir == "D":
            if self.mb >= self.mb_d:
                self.ballast_rate = 0

    def initialization(self):
        var = utils.load_json("vars/2d_glider_variables.json")
        pid_var = utils.load_json("vars/pid_variables.json")

        self.pid_control = var["pid_control"]
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

        self.theta_prev = pid_var["theta_prev"]

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

    def transformation(self):
        n2 = np.array([[self.phi, self.theta, self.psi]]).transpose()
        eta = np.array(
            [self.n1.transpose(), n2.transpose()], dtype=np.float32
        ).transpose()

        self.nu = np.array([self.v, self.Omega], dtype=np.float32)

        tau1 = self.F_ext
        tau2 = self.T_ext
        tau = np.array(
            [tau1.transpose(), tau2.transpose()], dtype=np.float32
        ).transpose()

        J1_n2, J2_n2 = utils.transformationMatrix(self.phi, self.theta, self.psi)

        self.n1_dot = J1_n2 * self.v
        self.n2_dot = J2_n2 * self.Omega

        R = J1_n2
        R_T = R.T  # same as np.linalg.inv(self.R)

        return R, R_T

    def control_transformation(self):

        if self.glide_dir == "D":
            self.wp = np.array(
                [[self.w1, self.controls.wp2, self.controls.wp3]]
            ).transpose()
        elif self.glide_dir == "U":
            self.wp = np.array(
                [[-self.w1, -self.controls.wp2, -self.controls.wp3]]
            ).transpose()

        self.wb = np.array([[0, 0, 0]]).transpose()

        self.ww = np.array([[0, 0, 0]]).transpose()

        # 3x3 matrix when mw != 0
        # self.F = np.array(
        #     [
        #         [
        #             np.linalg.inv(self.M)
        #             - self.rp_c @ (np.linalg.inv(self.J)) @ (self.rp_c)
        #             + (1 / self.mm) * self.I3,
        #             np.linalg.inv(self.M)
        #             - self.rp_c @ (np.linalg.inv(self.J)) @ (self.rb_c),
        #             np.linalg.inv(self.M)
        #             - self.rp_c @ (np.linalg.inv(self.J)) @ (self.rw_c)
        #         ],
        #         [
        #             np.linalg.inv(self.M)
        #             - self.rb_c @ (np.linalg.inv(self.J)) @ (self.rp_c),
        #             np.linalg.inv(self.M)
        #             - self.rb_c @ (np.linalg.inv(self.J)) @ (self.rb_c)
        #             + (1 / self.mb) * self.I3,
        #             np.linalg.inv(self.M)
        #             - self.rb_c @ (np.linalg.inv(self.J)) @ (self.rw_c)
        #         ],
        #         [
        #             np.linalg.inv(self.M)
        #             - self.rw_c @ (np.linalg.inv(self.J)) @ (self.rp_c),
        #             np.linalg.inv(self.M)
        #             - self.rw_c @ (np.linalg.inv(self.J)) @ (self.rb_c),
        #             np.linalg.inv(self.M)
        #             - self.rw_c @ (np.linalg.inv(self.J)) @ (self.rw_c)
        #             + (1 / self.mw) * self.I3,
        #         ]
        #     ]
        # )

        # 2x2 matrix when mw = 0
        self.F = np.array(
            [
                [
                    np.linalg.inv(self.M)
                    - self.rp_c @ (np.linalg.inv(self.J) @ self.rp_c)
                    + (1 / self.mm) * self.I3,
                    np.linalg.inv(self.M)
                    - self.rp_c @ (np.linalg.inv(self.J) @ self.rb_c),
                ],
                [
                    np.linalg.inv(self.M)
                    - self.rb_c @ (np.linalg.inv(self.J) @ self.rp_c),
                    np.linalg.inv(self.M)
                    - self.rb_c @ (np.linalg.inv(self.J) @ self.rb_c)
                    + (1 / self.mb) * self.I3,
                ],
            ]
        )

        self.H = np.linalg.inv(self.F)

        Zp = (
            -np.linalg.inv(self.M)
            @ (
                np.cross(self.M @ self.v + self.Pp + self.Pb, self.Omega, axis=0)
                + self.m0 * self.g * (self.R_T @ self.k_hat)
                + self.F_ext
            )
            - np.cross(self.Omega, self.rp_dot, axis=0)
            - np.linalg.inv(self.J)
            @ np.cross(
                np.cross(
                    self.J @ self.Omega + self.rp_c @ self.Pp + self.rb_c @ self.Pb,
                    self.Omega,
                    axis=0,
                )
                + np.cross(self.M @ self.v, self.v, axis=0)
                + self.T_ext
                + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
                + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
                + (self.mm * self.rp_c + self.mb * self.rb_c)
                * self.g
                @ (self.R_T @ self.k_hat),
                self.rp,
                axis=0,
            )
        )

        Zb = (
            -np.linalg.inv(self.M)
            @ (
                np.cross(self.M @ self.v + self.Pp + self.Pb, self.Omega, axis=0)
                + self.m0 * self.g * (self.R_T @ self.k_hat)
                + self.F_ext
            )
            - np.cross(self.Omega, self.rb_dot, axis=0)
            - np.linalg.inv(self.J)
            @ np.cross(
                np.cross(
                    self.J @ self.Omega + self.rp_c @ self.Pp + self.rb_c @ self.Pb,
                    self.Omega,
                    axis=0,
                )
                + np.cross(self.M @ self.v, self.v, axis=0)
                + self.T_ext
                + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
                + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
                + (self.mm * self.rp_c + self.mb * self.rb_c)
                * self.g
                @ (self.R_T @ self.k_hat),
                self.rb,
                axis=0,
            )
        )

        # Uncomment below if mw is not equal to 0
        # Zw = (
        #     -np.linalg.inv(self.M)
        #     @ (
        #         np.cross(self.M @ (self.v) + self.Pp + self.Pb, self.Omega, axis=0)
        #         + self.m0 * self.g * np.matmul(self.R_T, self.k_hat)
        #         + self.F_ext
        #     )
        #     - np.cross(self.Omega, self.rp_dot, axis=0)
        #     - np.cross(
        #         np.linalg.inv(self.J)
        #         @ (
        #             np.cross(
        #                 np.matmul(self.J, self.Omega)
        #                 + np.matmul(self.rp_c, self.Pp)
        #                 + np.matmul(self.rb_c, self.Pb),
        #                 self.Omega,
        #                 axis=0,
        #             )
        #             + np.cross(np.matmul(self.M, self.v), self.v, axis=0)
        #             + self.T_ext
        #             + np.cross(np.cross(self.Omega, self.rp, axis=0), self.Pp, axis=0)
        #             + np.cross(np.cross(self.Omega, self.rb, axis=0), self.Pb, axis=0)
        #             + (self.mm * self.rp_c + self.mb * self.rb_c)
        #             * self.g
        #             @ (self.R_T @ self.k_hat)
        #         ),
        #         self.rw,
        #         axis=0,
        #     )
        # )

        self.H = np.concatenate(
            (
                np.concatenate((self.H[0][0], self.H[0][1]), axis=0),
                np.concatenate((self.H[1][0], self.H[1][1]), axis=0),
            ),
            axis=1,
        )

        # self.u = self.H @ np.array(
        #     [(-Zp + self.wp).flatten(), (-Zb + self.wb).flatten(), (-Zw + self.ww).flatten()]
        # ).reshape(9,1) # uncomment if mw is not 0

        self.u = self.H @ np.array(((-Zp + self.wp), (-Zb + self.wb))).ravel()
        self.u = self.u.reshape(2, 3)

        self.u_bar = np.array([self.u[0]]).T

        self.u_b = np.array([self.Z3]).T

    def set_controls(self):
        self.control_transformation()

        if self.glide_dir == "D":
            if self.rp[0] >= self.rp1_d:
                self.u_bar = np.array([[0.0, 0.0, 0.0]]).T

        if self.glide_dir == "U":
            if self.rp[0] <= self.rp1_d:
                self.u_bar = np.array([[0.0, 0.0, 0.0]]).T

    def set_eom(self):
        if self.glide_dir == "U":
            self.Pp = self.mm * (
                self.v + np.cross(self.Omega, self.rp, axis=0) + self.rp_dot
            )

        elif self.glide_dir == "D":
            self.Pp = self.mm * (
                -self.v + np.cross(self.Omega, self.rp, axis=0) + self.rp_dot
            )

        self.Pb = np.array([[0.0, 0.0, 0.0]]).T

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

        n1_dot = self.R @ self.v

        Omega_dot = np.linalg.inv(self.J) @ T_bar

        v1_dot = np.linalg.inv(self.M) @ F_bar

        # Pp_dot = self.u_bar

        # Pb_dot = self.u_b

        # Pw_dot = self.u_w

        rp_ddot = self.wp

        rb_ddot = self.wb

        mb_dot = np.array([[self.ballast_rate, 0, 0]]).T

        return np.concatenate(
            [
                n1_dot,
                Omega_dot,
                v1_dot,
                self.rp_dot,
                self.rb_dot,
                # Pp_dot,
                # Pb_dot,
                rp_ddot,
                rb_ddot,
                mb_dot,
                np.array([self.n2_dot[1]]).T,
                self.wp,
            ]
        ).ravel()


if __name__ == "__main__":
    equations = Dynamics()
    equations.set_eom()
