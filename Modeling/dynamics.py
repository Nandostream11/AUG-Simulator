import numpy as np
import math
import utils
from glider_model import Vertical_Motion
from Parameters.slocum import SLOCUM_PARAMS

f_gravity = mt * g * np.dot(J1_n2.transpose(), k_hat)
f_buoyancy = -m * g * np.dot(J1_n2.transpose(), k_hat)
tau_gravity = np.cross(r_cg, mt * g * np.dot(J1_n2.transpose, k_hat))
tau_buoyancy = 0

# Terms to be defined:
# mb_d
# rp_d
# ballast_rate
# phi, theta, psi
# masses


class Dynamics:
    def __init__(self, z):
        # b same as n1, v same as v1, ohmega same as v2

        self.n1 = np.array(z[0:3])
        self.Omega = np.array(z[3:6])
        self.v1 = np.array(z[6:9])
        self.rp = np.array(z[9:12])
        self.rb = np.array(z[12:15])
        self.Pp = np.array(z[15:18])
        self.Pb = np.array(z[18:21])
        self.Pw = np.array(z[21:24])
        self.mb = z[24]
        self.theta = z[25]
        self.glide_dir = z[26]

        self.rp_c = utils.unit_vecs(self.rp)
        self.rb_c = utils.unit_vecs(self.rb)
        self.rw_c = np.array([0, 0, 0]).transpose()

        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        self.controls = SLOCUM_PARAMS.CONTROLS

        self.set_limits()

        self.transformation()

        self.set_force_torque()

    def transformation(self):
        n2 = np.array([[phi, self.theta, psi]]).transpose()
        eta = np.array(
            [self.n1.transpose(), n2.transpose()], dtype=np.float32
        ).transpose()
        self.v = np.array(
            [self.v1.transpose(), self.Omega.transpose()], dtype=np.float32
        ).transpose()
        # tau1 = np.array([[X, Y, Z]]).transpose()
        # tau2 = np.array([[K, M, N]]).transpose()
        # tau = np.array([tau1.transpose(), tau2.transpose()],
        #                dtype=np.float32).transpose()

        J1_n2, J2_n2 = utils.transformationMatrix(phi, self.theta, psi)

        # n1_dot = J1_n2 * self.v1
        # n2_dot = J2_n2 * self.Omega

        # self.eta_dot = np.array([n1_dot, n2_dot]).transpose()

        # self.eta_dot_T = np.linalg.inv(self.R)

        # Replace R with J1_n2 later

        self.R = J1_n2
        self.R_T = np.linalg.inv(self.R)

    def set_force_torque(self):
        self.alpha = np.arctan(self.v1[2] / self.v1[0])

        L = (KL0 + KL * self.alpha) * (
            math.pow(self.v1[0], 2) + math.pow(self.v1[2], 2)
        )
        D = (KD0 + KD * (math.pow(self.alpha, 2))) * (
            math.pow(self.v1[0], 2) + math.pow(self.v1[2], 2)
        )
        MDL = (KM0 + KM * self.alpha) * (
            math.pow(self.v1[0], 2)
            + math.pow(self.v1[2], 2)
            + KOmega1 * self.Omega[1]
            + KOmega2 * math.pow(self.Omega[1], 2)
        )

        self.F_ext = np.array([-D, 0, -L]).transpose()  # same as X, Y, Z

        self.T_ext = np.array([0, MDL, 0]).transpose()  # same as K, M, N

        return self.F_ext, self.T_ext

    def set_limits(self):
        if self.glide_dir == "D":
            if self.rp[0] >= rp_d:
                self.rp[0] = rp_d

            if self.mb >= self.mb_d:
                self.mb = self.mb_d
                ballast_rate = 0

            # Write for mw also

            self.rp_dot = (1 / mm) * self.Pp + self.v - np.cross(self.Omega, self.rp)

            # self.rb_dot = (1/mb) * self.Pb + self.v - np.cross(self.Omega, self.rb) # Uncomment if ballast mass moves
            # Assume ballast mass does not move
            self.rb_dot = self.Z3

            self.rw_dot = self.Z3

        elif self.glide_dir == "U":
            if self.rp[0] <= rp1_d:
                self.rp[0] = rp1_d

            if self.mb <= mb_d:
                self.mb = mb_d
                ballast_rate = 0

            # Write for mw also

            self.rp_dot = (1 / mm) * self.Pp - self.v - np.cross(self.Omega, self.rp)

            # self.rb_dot = (1/mb) * self.Pb - self.v - np.cross(self.Omega, self.rb) # Uncomment if ballast mass moves
            self.rb_dot = self.Z3

            self.rw_dot = self.Z3

    def control_transformation(self):
        if self.glide_dir == "D":
            self.wp = np.array(
                [self.controls.wp1, self.controls.wp2, self.controls.wp3]
            ).transpose()
        elif self.glide_dir == "U":
            self.wp = np.array(
                [-self.controls.wp1, self.controls.wp2, self.controls.wp3]
            ).transpose()

        self.wb = np.array(
            [self.controls.wb1, self.controls.wb2, self.controls.wb3]
        ).transpose()

        self.ww = np.array(
            [self.controls.ww1, self.controls.ww2, self.controls.ww3]
        ).transpose()

        self.F = np.array(
            [
                [
                    np.linalg.inv(M)
                    - self.rp_c * np.linalg.inv(J) * self.rp_c
                    + (1 / mm) * self.I3,
                    np.linalg.inv(M) - self.rp_c * np.linalg.inv(J) * self.rb_c,
                    np.linalg.inv(M) - self.rp_c * np.linalg.inv(J) * self.rw_c,
                ],
                [
                    np.linalg.inv(M) - self.rp_c * np.linalg.inv(J) * self.rp_c,
                    np.linalg.inv(M)
                    - self.rp_c * np.linalg.inv(J) * self.rb_c
                    + (1 / mb) * self.I3,
                    np.linalg.inv(M) - self.rp_c * np.linalg.inv(J) * self.rw_c,
                ],
                [
                    np.linalg.inv(M) - self.rp_c * np.linalg.inv(J) * self.rp_c,
                    np.linalg.inv(M) - self.rp_c * np.linalg.inv(J) * self.rb_c,
                    np.linalg.inv(M)
                    - self.rp_c * np.linalg.inv(J) * self.rw_c
                    + (1 / mw) * self.I3,
                ],
            ]
        )

        self.H = np.linalg.inv(self.F)

        Zp = (
            -np.linalg.inv(M)
            * (
                np.cross(M * self.v1 + self.Pp + self.Pb + self.Pw, self.Omega)
                + m0 * self.g * self.R_T * self.k_hat
                + self.F_ext
            )
            - np.cross(self.Omega, self.rp_dot)
            - np.linalg.inv(J)
            * (
                np.cross(
                    J * self.Omega
                    + self.rp_c * self.Pp
                    + self.rb_c * self.Pb
                    + self.rw_c * self.Pw,
                    self.Omega,
                )
                + np.cross(M * self.v1, self.v1)
                + self.T_ext
            )
            + np.cross(np.cross(self.Omega, self.rp), self.Pp)
            + np.cross(np.cross(self.Omega, self.rb), self.Pb)
            + np.cross(np.cross(self.Omega, self.rw), self.Pw)
            + np.cross(
                (mm * self.rp_c + mb * self.rb_c + mw * self.rw_c)
                * self.R_T
                * self.k_hat,
                self.rp,
            )
        )

        Zb = (
            -np.linalg.inv(M)
            * (
                np.cross(M * self.v1 + self.Pp + self.Pb + self.Pw, self.Omega)
                + m0 * self.g * self.R_T * self.k_hat
                + self.F_ext
            )
            - np.cross(self.Omega, self.rb_dot)
            - np.linalg.inv(J)
            * (
                np.cross(
                    J * self.Omega
                    + self.rp_c * self.Pp
                    + self.rb_c * self.Pb
                    + self.rw_c * self.Pw,
                    self.Omega,
                )
                + np.cross(M * self.v1, self.v1)
                + self.T_ext
            )
            + np.cross(np.cross(self.Omega, self.rp), self.Pp)
            + np.cross(np.cross(self.Omega, self.rb), self.Pb)
            + np.cross(np.cross(self.Omega, self.rw), self.Pw)
            + np.cross(
                (mm * self.rp_c + mb * self.rb_c + mw * self.rw_c)
                * self.R_T
                * self.k_hat,
                self.rb,
            )
        )

        Zw = (
            -np.linalg.inv(M)
            * (
                np.cross(M * self.v1 + self.Pp + self.Pb + self.Pw, self.Omega)
                + m0 * self.g * self.R_T * self.k_hat
                + self.F_ext
            )
            - np.cross(self.Omega, self.rw_dot)
            - np.linalg.inv(J)
            * (
                np.cross(
                    J * self.Omega
                    + self.rp_c * self.Pp
                    + self.rb_c * self.Pb
                    + self.rw_c * self.Pw,
                    self.Omega,
                )
                + np.cross(M * self.v1, self.v1)
                + self.T_ext
            )
            + np.cross(np.cross(self.Omega, self.rp), self.Pp)
            + np.cross(np.cross(self.Omega, self.rb), self.Pb)
            + np.cross(np.cross(self.Omega, self.rw), self.Pw)
            + np.cross(
                (mm * self.rp_c + mb * self.rb_c + mw * self.rw_c)
                * self.R_T
                * self.k_hat,
                self.rw,
            )
        )

        self.u = self.H * np.array([[-Zp + self.wp], [-Zb + self.wb], [-Zw + self.ww]])

        self.u_b = self.Z3
        self.u_w = self.Z3

    def set_controls(self):
        self.control_transformation()

        if self.glide_dir == "D":
            if self.rp[0] >= rp1_d:
                self.u_bar = self.Z3
            else:
                self.u_bar = self.u[0]

        if self.glide_dir == "U":
            if self.rp[0] <= rp1_d:
                self.u_bar = self.Z3
            else:
                self.u_bar = self.u[0]

        # Write for u_b and u_w as well

    def set_eom(self):
        self.set_controls()

        T_bar = (
            np.cross(
                (
                    J * self.Omega
                    + self.rp_c * self.Pp
                    + self.rb_c * self.Pb
                    + self.rw_c * self.Pw
                ),
                self.Omega,
            )
            + np.cross(M * self.v1, self.v1)
            + np.cross(np.cross(self.Omega, self.rp), self.Pp)
            + np.cross(np.cross(self.Omega, self.rb), self.Pb)
            + np.cross(np.cross(self.Omega, self.rw), self.Pw)
            + (mm * self.rp_c + mb * self.rb_c + mw * self.rw_c)
            * self.g
            * self.R_T
            * self.k_hat
            + self.T_ext
            - self.rp_c * self.u_bar
            - (self.rb_c * self.u_b + self.rw_c * self.u_w)
        )

        F_bar = (
            np.cross((M * self.v1 + self.Pp + self.Pb + self.Pw), self.Omega)
            + m0 * self.g * self.R_T * self.k_hat
            + self.F_ext
            - self.u_bar
            - (self.u_b + self.u_w)
        )

        # R_dot

        n1_dot = self.R * self.v1

        Omega_dot = np.linalg.inv(J) * T_bar

        v1_dot = np.linalg.inv(M) * F_bar

        # rp_dot solved earlier

        # rb_dot solved earlier

        # rw_dot solved earlier

        Pp_dot = self.u_bar

        Pb_dot = self.u_b

        Pw_dot = self.u_w

        mb_dot = ballast_rate

        theta_dot = self.Omega[1]

        return [
            n1_dot.transpose(),
            Omega_dot.transpose(),
            v1_dot.transpose(),
            self.rp_dot.transpose(),
            self.rb_dot.transpose(),
            self.rw_dot.transpose(),
            Pp_dot.transpose(),
            Pb_dot.transpose(),
            Pw_dot.transpose,
            mb_dot,
            theta_dot,
        ]


if __name__ == "__main__":
    equations = Dynamics()
    equations.set_eom()
