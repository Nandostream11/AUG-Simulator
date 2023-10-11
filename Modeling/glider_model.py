import numpy as np
import math
from scipy.integrate import solve_ivp
from scipy.integrate import odeint
from scipy.integrate import RK45
from scipy.integrate import ode
from scipy import integrate
import utils
from Parameters.slocum import SLOCUM_PARAMS
from Modeling.dynamics import Dynamics
import matplotlib.pyplot as plt


class Vertical_Motion:
    def __init__(self):
        self.mass_params = SLOCUM_PARAMS.GLIDER_CONFIG
        self.hydro_params = SLOCUM_PARAMS.HYDRODYNAMICS
        self.vars = SLOCUM_PARAMS.VARIABLES

        self.initialization("SLOCUM")

        self.solver_array = []
        self.total_time = []
        self.derivative = []

    def initialization(self, glider_name="SLOCUM"):
        self.g, self.I3, self.Z3, self.i_hat, self.j_hat, self.k_hat = utils.constants()

        if glider_name == "SLOCUM":
            from Parameters.slocum import SLOCUM_PARAMS as P

        self.mass_params = P.GLIDER_CONFIG
        self.hydro_params = P.HYDRODYNAMICS
        self.vars = P.VARIABLES

        self.mh = 40.0
        self.mw = 0.0
        self.mb = 1.0
        self.mm = 9.0

        self.ms = 0

        self.mt = self.mh + self.mw + self.mb + self.ms + self.mm

        # self.ms = self.mh + self.mw + self.mb
        # self.mt = self.ms + self.mm

        self.m = 50.0
        self.m0 = self.mt - self.m

        self.Mf = np.diag(
            [5, 60, 70]
        )
        self.Jf = np.diag(
            [4, 12, 11]
        )

        self.M = self.mh * self.I3 + self.Mf
        self.J = self.Jf  # J = Jf + Jh

        self.KL = 132.5
        self.KL0 = 0.0
        self.KD = 25
        self.KD0 = 2.15
        self.KM = -100
        self.KM0 = 0.0
        self.KOmega1 = -50
        self.KOmega2 = -50
        
        self.rp3 = 0.05
        self.rs3 = 0.0
        self.rb1 = 0.0
        self.rb3 = 0.0

        # [self.rp1, self.rp2, self.rp3] = [0.0, 0.0, 0.05]
        # [self.rb1, self.rb2, self.rb3] = [0.0, 0.0, 0.0]
        # [self.rw1, self.rw2, self.rw3] = [0.0, 0.0, 0.0]
        # [self.rs1, self.rs2, self.rs3] = [0.0, 0.0, 0.0]


        self.glide_angle_deg = 25
        self.V_d = 0.3
        self.ballast_rate = 0.025

        self.set_first_run_params()

    def set_first_run_params(self):
        # self.n1_0 = np.array([[0.0, 0.0, 0.0]])
        # self.Omega0 = np.array([[0.0, 0.0, 0.0]])
        self.phi = 0
        self.theta0 = math.radians(25)
        self.psi = 0

        # self.Pp = [0.0, 0.0, 0.0]
        # # self.Pb = [0.0, 0.0, 0.0]
        # self.Pw = [0.0, 0.0, 0.0]
        # self.Ps = [0.0, 0.0, 0.0]

    def set_desired_trajectory(self):
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
        for i in range(l):
            self.e_i_d = self.E_i_d[i]

            print(
                "Iteration {} | Desired glide angle in deg = {}".format(
                    i, math.degrees(self.e_i_d)
                )
            )

            if (self.e_i_d) > 0:
                self.glider_direction = "U"
                self.ballast_rate = -abs(self.ballast_rate)
                print("Glider moving in upward direction\n")

            elif (self.e_i_d) < 0:
                self.glider_direction = "D"
                self.ballast_rate = abs(self.ballast_rate)
                print("Glider moving in downward direction\n")

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

            self.mb_d = (self.m - self.mh - self.mm - self.mw - self.ms) + (
                1 / self.g
            ) * (
                -math.sin(self.e_i_d) * (self.KD0 + self.KD * math.pow(self.alpha_d, 2))
                + math.cos(self.e_i_d) * (self.KL0 + self.KL * self.alpha_d)
            ) * math.pow(
                self.V_d, 2
            )

            self.theta_d = self.e_i_d + self.alpha_d
            lambda_val = 1

            self.v1_d = self.V_d * math.cos(self.alpha_d)
            self.v3_d = self.V_d * math.sin(self.alpha_d)

            # self.Pp1_d = self.mm * self.v1_d
            # self.Pp3_d = self.mm * self.v3_d

            # self.Pb1_d = self.mb * self.v1_d
            # self.Pb3_d = self.mb * self.v3_d

            # self.m0_d = self.mb_d + self.mh + self.mm - self.m

            # self.rp1_d = -self.rp3 * math.tan(self.theta_d) + (
            #     1 / (self.mm * self.g * math.cos(self.theta_d))
            # ) * (
            #     (self.Mf[2, 2] - self.Mf[0, 0]) * self.v1_d * self.v3_d
            #     + (self.KM0 + self.KM * self.alpha_d) * math.pow(self.V_d, 2)
            # )

            self.rp1_d = (1 / (self.mm + lambda_val * self.ms)) * (
                -self.mb_d * self.rb1
                - math.tan(self.theta_d)
                * (self.mm * self.rp3 + self.ms * self.rs3 + self.mb * self.rb3)
                + (1 / (self.g * math.cos(self.theta_d)))
                * (
                    (self.Mf[2, 2] - self.Mf[0, 0]) * self.v1_d * self.v3_d
                    + (self.KM0 + self.KM * self.alpha_d) * math.pow(self.V_d, 2)
                )
            )
            self.rs1_d = lambda_val * self.rp1_d
            
            self.save_json()

            # These are the initial conditions at every peak of the sawtooth trajectory
            
            
            if i == 0:
                self.z_in = np.array([[0.0, 0.0, 0.0],
                                      [0.0, 0.0, 0.0],
                                      [self.v1_d, 0.0, self.v3_d],
                                      [self.rp1_d, 0.0, 0.05],
                                      [self.rs1_d, 0.0, self.rs3],
                                      [self.rb1, 0.0, self.rb3],
                                      [0.0, 0.0, 0.0],
                                      [0.0, 0.0, 0.0],
                                      [self.mb * self.v1_d, 0.0, self.mb * self.v3_d],
                                      [self.mb_d, 0, 0],
                                      [self.theta0, 0, 0]]).flatten()
            
            else:
                p = len(self.solver_array)
                self.z_in = np.array([[self.solver_array[-1][0:3]],
                                      [self.solver_array[-1][3:6]],
                                      [self.solver_array[-1][6:9]],
                                      [self.solver_array[-1][9:12]],
                                      [self.solver_array[-1][12:15]],
                                      [self.solver_array[-1][15:18]],
                                      [self.solver_array[-1][18:21]],
                                      [self.solver_array[-1][21:24]],
                                      [self.solver_array[-1][24:27]],
                                      [self.solver_array[-1][27:30]],
                                      [self.solver_array[-1][30:33]]]).flatten()
            

            self.i = i

            sol = self.solve_ode(self.z_in)  
            # sol.y[11] = 0.05    
            
            if i == 0:
                self.solver_array = np.array(sol.y).T
                self.total_time = np.array(sol.t)
            else:
                self.solver_array = np.concatenate((self.solver_array, np.array(sol.y).T))
                self.total_time = np.concatenate((self.total_time, np.array(sol.t)))
                
            if i == l - 1:
                self.plots()                
                

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
            # "m0_d": self.m0_d,
            "rp1_d": self.rp1_d,
            "rs1_d": self.rs1_d,
            "rp3": 0.05,
            "rb1": self.rb1,
            "rb3": self.rb3,
            # "rw1": self.rw1,
            # "rw2": self.rw2,
            # "rw3": self.rw3,
            "rs3": self.rs3,
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
        glide_vars["ms"] = self.ms
        glide_vars["m"] = self.m
        glide_vars["m0"] = self.m0
        glide_vars["mt"] = self.mt

        utils.save_json(glide_vars)

    def solve_ode(self, z0):
        
        # eom = Dynamics(z0)
        
        def dvdt(t, y):
            eom = Dynamics(y)
            D = eom.set_eom()
            self.derivative = np.concatenate((self.derivative, np.array([D[11]])))
            return D

        self.t = np.linspace(375 * (self.i), 375 * (self.i + 1))[:-1]
        # self.t = np.arange(200 * self.i, 200 * (self.i + 1))
        
        # sol = solve_ivp(
        #     dvdt, t_span=(0, max(self.t)), y0=z0, method='RK45', t_eval=self.t, dense_output=True)

        sol = solve_ivp(
            dvdt, t_span=(min(self.t), max(self.t)), y0=z0, method='RK45', t_eval=self.t, dense_output=True)
        
        # sol = odeint(dvdt, y0=z0, t=self.t, tfirst=True)
        # y = sol
        # t = self.t
        
        # plt.plot(sol.t, sol.y[2])
        # plt.show()
        
        return sol
    
    def new_ode_func(self):
        
        def dvdt(t, y):
            eom = Dynamics(y)
            D = eom.set_eom()
            Z = concat(np.array(D, dtype="object"))
            return Z

        def concat(z_in):
            z_in = np.concatenate(
                (
                    z_in[0].ravel(),
                    z_in[1].ravel(),
                    z_in[2].ravel(),
                    z_in[3].ravel(),
                    z_in[4].ravel(),
                    z_in[5].ravel(),
                    z_in[6].ravel(),
                    z_in[7].ravel(),
                    z_in[8].ravel(),
                    z_in[9].ravel(),
                    z_in[10].ravel(),
                )
            )
            return z_in
        
        self.z0 = concat(self.z_in)
        y = self.z0
        
        t = np.linspace(375 * (self.i), 375 * (self.i + 1))
        t0 = t[0]
        tf = t[-1]
        meth = 'bdf'
        time = [t0, tf]
        r = integrate.solve_ivp(fun=lambda t,y: dvdt(t,y), t_span=time, y0=y, method='RK45', dense_output=True, rtol=1e-13, atol=1e-22)
        k = 0
        while r.t[k] < tf:
            r = integrate.solve_ivp(fun=lambda t,y: dvdt(t,y), t_span=time, y0=y, method='RK45', dense_output=True, rtol=1e-13, atol=1e-22)
            print(r.t[k])
            k += 1
        new_time = r.t
        new_temp = r.y[:, -1]
        
        r = ode(dvdt)
        r.set_integrator('vode', method='bdf', rtol=1e-13, atol=1e-22, with_jacobian=False)
        r.set_initial_value(y, t0)
        r.integrate(tf)
        
        t = []
        Y = [y]
        
        while r.t < tf:
            r.integrate(tf, step=True)
            Y = np.vstack((Y, [r.y]))
            t.append([r.t])
            
        new_temp = Y[-1, :]
        

    def plots(self):
        # breakpoint()
        plt.plot(self.total_time, self.solver_array.T[11])  
        # plt.plot(np.linspace(0,5312, 5312), self.derivative)      
        plt.show()


if __name__ == "__main__":
    Z = Vertical_Motion()
    Z.set_desired_trajectory()
    # Z.save_json()
