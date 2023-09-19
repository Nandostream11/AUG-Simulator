import numpy as np
import math
import utils
import os
import pdb
from Parameters.slocum import SLOCUM_PARAMS


class Vertical_Motion:
    def __init__(self):            

        self.mass_params = SLOCUM_PARAMS.GLIDER_CONFIG
        self.hydro_params = SLOCUM_PARAMS.HYDRODYNAMICS
        self.vars = SLOCUM_PARAMS.VARIABLES
                
        self.initialization()
                
    def constants(self):
        self.g = 9.816
        self.I3 = np.eye(3)
        self.Z3 = np.zeros(3)

    def initialization(self):
        
        self.constants()

        self.mh = self.mass_params.HULL_MASS
        self.mw = self.mass_params.FIXED_POINT_MASS
        self.mb = self.mass_params.BALLAST_MASS
        self.mm = self.mass_params.INT_MOVABLE_MASS

        self.ms = self.mh + self.mw + self.mb
        self.mt = self.ms + self.mm

        self.m = self.mass_params.FLUID_DISP_MASS
        self.m0 = self.mt - self.m

        self.Mf = np.diag(
            [self.mass_params.MF1, self.mass_params.MF2, self.mass_params.MF3])
        self.Jf = np.diag(
            [self.mass_params.J1, self.mass_params.J2, self.mass_params.J3])

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

        [self.rp1, self.rp2, self.rp3] = [
            self.vars.rp1, self.vars.rp2, self.vars.rp3]
        [self.rb1, self.rb2, self.rb3] = [
            self.vars.rb1, self.vars.rb2, self.vars.rb3]

        self.glide_angle_deg = self.vars.GLIDE_ANGLE
        self.V_d = self.vars.SPEED
        self.u_ballast_rate = self.vars.BALLAST_RATE
        
        self.set_first_run_params()
        
    def set_first_run_params(self):
        self.b0 = np.array([0.0, 0.0, 0.0]).transpose()
        self.Omega0 = np.array([0.0, 0.0, 0.0]).transpose()
        self.phi0 = self.mass_params.PHI
        self.theta0 = self.mass_params.THETA
        self.psi0 = self.mass_params.PSI
        
        self.Pp = np.array([0.0, 0.0, 0.0]).transpose()
        self.Pb = np.array([0.0, 0.0, 0.0]).transpose()
        self.Pw = np.array([0.0, 0.0, 0.0]).transpose()
        

    def set_desired_trajectory(self):

        self.E_i_d = np.array([math.radians(-self.glide_angle_deg), math.radians(self.glide_angle_deg),
                               math.radians(-self.glide_angle_deg), math.radians(self.glide_angle_deg)])

        lim1 = math.degrees(math.atan(2*(self.KD/self.KL)*(self.KL0/self.KL +
                            math.sqrt(math.pow(self.KL0/self.KL, 2) + self.KD0/self.KD))))

        lim2 = math.degrees(math.atan(2*(self.KD/self.KL)*(self.KL0/self.KL -
                            math.sqrt(math.pow(self.KL0/self.KL, 2) + self.KD0/self.KD))))

        l = len(self.E_i_d)

        for i in range(l):

            e_i_d = self.E_i_d[i]

            print('Iteration {} | Desired glide angle in deg = {}'.format(i, math.degrees(e_i_d)))

            if e_i_d > 0:
                glider_direction = 'U'
                u_ballast_rate = -abs(self.u_ballast_rate)
                print('Glider moving in upward direction\n')

            elif e_i_d < 0:
                glider_direction = 'D'
                u_ballast_rate = abs(self.u_ballast_rate)
                print('Glider moving in downward direction\n')

            alpha_d = (1/2)*(self.KL/self.KD)*math.tan(e_i_d)*(-1 + math.sqrt(1 - 4*(self.KD /
                                                                                     math.pow(self.KL, 2))*(1/math.tan(e_i_d))*(self.KD0*(1/math.tan(e_i_d)) + self.KL0)))

            mb_d = (self.m - self.mh - self.mm) + (1/self.g) * (-math.sin(e_i_d)*(self.KD0 + self.KD *
                                                                                  math.pow(alpha_d, 2)) + math.cos(e_i_d)*(self.KL0 + self.KL*alpha_d)) * math.pow(self.V_d, 2)

            theta_d = e_i_d + alpha_d

            v1_d = self.V_d * math.cos(alpha_d)
            v3_d = self.V_d * math.sin(alpha_d)

            Pp1_d = self.mm * v1_d
            Pp3_d = self.mm * v3_d
            
            Pb1_d = self.mb * v1_d
            Pb3_d = self.mb * v3_d            

            m0_d = mb_d + self.mh + self.mm - self.m

            self.rp1_d = -self.rp3 * math.tan(theta_d) + (1/(self.mm*self.g*math.cos(theta_d))) * \
                ((self.Mf[2]-self.Mf[0])*v1_d*v3_d +
                 (self.KM0+self.KM*alpha_d)*math.pow(self.V_d, 2))
                
            
            if i == 0:
                
                Z = [self.b0,
                     self.Omega0,
                     [v1_d, 0.0, v3_d],
                     [self.rp1_d, self.rp2, self.rp3],
                     [self.rb1, self.rb2, self.rb3],
                     [Pp1_d, self.Pp[1], Pp3_d],
                     [Pb1_d, self.Pb[1], Pb3_d],
                     self.Pw,
                     mb_d,
                     self.theta0]
            
            else:
                
                Z = None                
        
    def plots(self):
        pass


if __name__ == '__main__':
    
    Z = Vertical_Motion()
    Z.set_desired_trajectory()
