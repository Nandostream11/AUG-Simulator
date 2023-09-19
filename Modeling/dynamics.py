import numpy as np
import math
import utils
from glider_model import Vertical_Motion

class Dynamics:
    def __init__(self):
        pass


# Physics
g = 9.816 # acceleration due to gravity

# Vectors
[x, y, z] = [0.0, 0.0, 0.0]
[phi, theta, psi] = [0.0, 0.0, 0.0]
[u, v, w] = [0.0, 0.0, 0.0]
[p, q, r] = [0.0, 0.0, 0.0]
[X, Y, Z] = [0.0, 0.0, 0.0]
[K, M, N] = [0.0, 0.0, 0.0]

i_hat = np.array([[1,0,0]]).transpose()
j_hat = np.array([[0,1,0]]).transpose()
k_hat = np.array([[0,0,1]]).transpose()


n1 = np.array([[x, y, z]]).transpose()
n2 = np.array([[phi, theta, psi]]).transpose()
eta = np.array([n1.transpose(), n2.transpose()], dtype=np.float32).transpose()
v1 = np.array([[u, v, w]]).transpose()
v2 = np.array([[p, q, r]]).transpose()
v = np.array([v1.transpose(), v2.transpose()], dtype=np.float32).transpose()
tau1 = np.array([[X, Y, Z]]).transpose()
tau2 = np.array([[K, M, N]]).transpose()
tau = np.array([tau1.transpose(), tau2.transpose()],
               dtype=np.float32).transpose()

J1_n2, J2_n2 = utils.transformationMatrix(phi, theta, psi)

n1_dot = J1_n2 * v1
n2_dot = J2_n2 * v2

eta_dot = np.array([n1_dot, n2_dot]).transpose()


# Vehicle Model
mh = 0.0  # Hull mass
mw = 0.0  # Fixed point mass
mb = 0.0  # variable ballast point mass
ms = mh + mw + mb  # stationary mass of glider
mm = 0.0  # internal movable mass

mt = ms + mm  # total mass of UG

rho_a = 0.0 # mass density of vehicle

r_h = 0.0  # zero since hull mass is uniformly distributed
r_w = 0.0
r_b = 0.0
# r_p = np.array([[rp1, rp2, rp3]]).transpose()

m = 0.0 # mass of fluid displaced by vehicle
m0 = mt - m # net buoyancy

r_cg = (mh*r_h + mw*r_w + mb*r_b + mm*r_p) / mt

f_gravity = mt * g * np.dot(J1_n2.transpose(),k_hat)
f_buoyancy = -m * g * np.dot(J1_n2.transpose(),k_hat)
tau_gravity = np.cross(r_cg, mt*g*np.dot(J1_n2.transpose,k_hat))
tau_buoyancy = 0 

# Jh = none # inertia matrix of uniformly distributed hull mass
# Js = Jh - mw * (r_w/) 

ed = 0 # desired glide path angle
Vd = 0 # desired speed of vehicle