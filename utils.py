import numpy as np
import math
import json


def transformationMatrix(phi, theta, psi):
    J1_n2 = np.array(
        [
            [
                math.cos(phi) * math.cos(theta),
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
    # mag = math.sqrt(math.pow(r[0], 2) + math.pow(r[1], 2) + math.pow(r[2], 2))

    # return (1 / mag) * r

    return np.array([[0, -r[2], r[1]], [r[2], 0, -r[0]], [-r[1], r[0], 0]])


def constants():
    g = 9.816
    I3 = np.eye(3)
    Z3 = np.zeros(3)

    i_hat = np.array([1, 0, 0]).transpose()
    j_hat = np.array([0, 1, 0]).transpose()
    k_hat = np.array([0, 0, 1]).transpose()

    return g, I3, Z3, i_hat, j_hat, k_hat


def save_json(vars, path="vars/glider_variables.json"):
    with open(path, "w", encoding="utf-8") as file:
        json.dump(vars, file, separators=(",", ":"), sort_keys=True, indent=4)


def load_json(path="vars/glider_variables.json"):
    file = open(path)

    return json.load(file)
