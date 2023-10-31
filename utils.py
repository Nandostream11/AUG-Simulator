import numpy as np
import math
import json
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def transformationMatrix(phi, theta, psi):
    J1_n2 = np.array(
        [
            [
                math.cos(psi) * math.cos(theta),
                -math.sin(psi) * math.cos(phi)
                + math.cos(psi) * math.sin(theta) * math.sin(phi),
                math.sin(psi) * math.sin(phi)
                + math.cos(psi) * math.cos(phi) * math.sin(theta),
            ],
            [
                math.sin(psi) * math.cos(theta),
                math.cos(psi) * math.cos(phi)
                + math.sin(phi) * math.sin(theta) * math.sin(psi),
                -math.cos(psi) * math.sin(phi)
                + math.sin(theta) * math.sin(psi) * math.cos(phi),
            ],
            [
                -math.sin(theta),
                math.cos(theta) * math.sin(phi),
                math.cos(theta) * math.cos(phi),
            ],
        ]
    )

    J2_n2 = np.array(
        [
            [1, math.sin(theta) * math.tan(theta), math.cos(phi) * math.tan(theta)],
            [0, math.cos(phi), -math.sin(phi)],
            [0, math.sin(phi) / math.cos(theta), math.cos(phi) / math.cos(theta)],
        ]
    )

    return J1_n2, J2_n2


def unit_vecs(r):
    return np.array(
        [[0, -r[2][0], r[1][0]], [r[2][0], 0, -r[0][0]], [-r[1][0], r[0][0], 0]]
    )  # Skew-symmetric matrix


def constants():
    g = 9.816
    I3 = np.eye(3)
    Z3 = np.zeros(3)

    i_hat = np.array([[1, 0, 0]]).transpose()
    j_hat = np.array([[0, 1, 0]]).transpose()
    k_hat = np.array([[0, 0, 1]]).transpose()

    return g, I3, Z3, i_hat, j_hat, k_hat


def save_json(vars, path="vars/2d_glider_variables.json"):
    with open(path, "w", encoding="utf-8") as file:
        json.dump(vars, file, separators=(",", ":"), sort_keys=True, indent=4)


def load_json(path="vars/2d_glider_variables.json"):
    file = open(path)

    return json.load(file)


def plots(t, x, plot):
    vel = []
    phi = []
    theta = []
    psi = []
    for i in range(len(t)):
        v = math.sqrt(math.pow(x[6][i], 2) + math.pow(x[8][i], 2))
        phi_angle = math.degrees(x[-3][i])
        theta_angle = math.degrees(x[-2][i])
        psi_angle = math.degrees(x[-1][i]) % 360
        vel.append(v)
        phi.append(phi_angle)
        theta.append(theta_angle)
        psi.append(psi_angle)

    vel = np.array(vel)
    phi = np.array(phi)
    theta = np.array(theta)
    psi = np.array(psi)

    if plot == ["3D"]:
        ax = plt.axes(projection="3d")
        ax.plot3D(x[0], x[1], x[2], "gray")
        ax.set_xlabel("x (m)")
        ax.set_ylabel("y (m)")
        ax.set_zlabel("z (m)")
        ax.invert_zaxis()
        plt.show()

    elif plot == ["all"] or plot == "all":
        fig, ax = plt.subplots(3, 2)
        ax[0, 0].plot(t, x[0])
        ax[0, 0].set(xlabel="time (s)", ylabel="x (m)")
        ax[1, 0].plot(t, x[1])
        ax[1, 0].set(xlabel="time (s)", ylabel="y (m)")
        ax[2, 0].plot(t, x[2])
        ax[2, 0].set(xlabel="time (s)", ylabel="z (m)")
        ax[0, 1].plot(t, phi)
        ax[0, 1].set(xlabel="time (s)", ylabel="phi (deg)")
        ax[1, 1].plot(t, theta)
        ax[1, 1].set(xlabel="time (s)", ylabel="theta (deg)")
        ax[2, 1].plot(t, psi)
        ax[2, 1].set(xlabel="time (s)", ylabel="psi (deg)")
        plt.show()

        fig, ax = plt.subplots(4, 2)
        ax[0, 0].plot(t, x[3])
        ax[0, 0].set(xlabel="time (s)", ylabel="Omega1 (rad/s)")
        ax[1, 0].plot(t, x[4])
        ax[1, 0].set(xlabel="time (s)", ylabel="Omega2 (rad/s)")
        ax[2, 0].plot(t, x[5])
        ax[2, 0].set(xlabel="time (s)", ylabel="Omega3 (rad/s)")
        ax[0, 1].plot(t, x[6])
        ax[0, 1].set(xlabel="time (s)", ylabel="v1 (m/s)")
        ax[1, 1].plot(t, x[7])
        ax[1, 1].set(xlabel="time (s)", ylabel="v2 (m/s)")
        ax[2, 1].plot(t, x[8])
        ax[2, 1].set(xlabel="time (s)", ylabel="v3 (m/s)")
        ax[3, 0].plot(t, vel)
        ax[3, 0].set(xlabel="time (s)", ylabel="velocity (m/s)")
        plt.show()

        fig, ax = plt.subplots(2, 2)
        ax[0, 0].plot(t, x[9])
        ax[0, 0].set(xlabel="time (s)", ylabel="rp1 (m)")
        ax[1, 0].plot(t, x[10])
        ax[1, 0].set(xlabel="time (s)", ylabel="rp2 (m)")
        ax[0, 1].plot(t, x[11])
        ax[0, 1].set(xlabel="time (s)", ylabel="rp3 (m)")
        ax[1, 1].plot(t, x[21])
        ax[1, 1].set(xlabel="time (s)", ylabel="mb (kg)")
        plt.show()

    else:
        for p in plot:
            if p == "x":
                plt.plot(t, x[0])
            elif p == "y":
                plt.plot(t, x[1])
            elif p == "z":
                plt.plot(t, x[2])
            elif p == "omega1":
                plt.plot(t, x[3])
            elif p == "omega2":
                plt.plot(t, x[4])
            elif p == "omega3":
                plt.plot(t, x[5])
            elif p == "v1":
                plt.plot(t, x[6])
            elif p == "v2":
                plt.plot(t, x[7])
            elif p == "v3":
                plt.plot(t, x[8])
            elif p == "vel":
                plt.plot(t, vel)
            elif p == "rp1":
                plt.plot(t, x[9])
            elif p == "rp2":
                plt.plot(t, x[10])
            elif p == "rp3":
                plt.plot(t, x[11])
            elif p == "mb":
                plt.plot(t, x[21])
            elif p == "phi":
                plt.plot(t, phi)
            elif p == "theta":
                plt.plot(t, theta)
            elif p == "psi":
                plt.plot(t, psi)

            plt.show()
