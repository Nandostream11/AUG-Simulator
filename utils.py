import numpy as np
import math


def transformationMatrix(phi, theta, psi):
    J1_n2 = np.array([[math.cos(phi)*math.cos(theta), -math.sin(psi)*math.cos(phi)+math.cos(psi)*math.sin(theta)
                       * math.sin(phi), math.sin(psi)*math.sin(phi)+math.cos(psi)*math.cos(phi)*math.sin(theta)],
                      [math.sin(psi)*math.cos(theta), math.cos(psi)*math.cos(phi)+math.sin(phi)*math.sin(theta)
                       * math.sin(psi), -math.cos(psi)*math.sin(phi)+math.sin(theta)*math.sin(psi)*math.cos(phi)],
                      [-math.sin(theta), math.cos(theta)*math.sin(phi), math.cos(theta)*math.cos(phi)]])

    J2_n2 = np.array([[1, math.sin(theta)*math.tan(theta), math.cos(phi)*math.tan(theta)],
                      [0, math.cos(phi), -math.sin(phi)],
                      [0, math.sin(phi)/math.cos(theta), math.cos(phi)/math.cos(theta)]])

    return J1_n2, J2_n2
