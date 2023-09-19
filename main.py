import os
import pdb
from Modeling.glider_model import Vertical_Motion

def main():
    Z = Vertical_Motion()
    Z.set_desired_trajectory()

if __name__ == '__main__':
    main()