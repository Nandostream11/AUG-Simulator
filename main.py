import os
import pdb
from Modeling.glider_model import Vertical_Motion
from Modeling.dynamics import Dynamics
import numpy as np


def main():
    Z = Vertical_Motion()
    Z.set_desired_trajectory()
    # x = np.zeros(30)
    
    # D = Dynamics(x)
    # D.initialization()


if __name__ == "__main__":
    main()
