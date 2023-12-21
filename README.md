# Autonomous Underwater Glider Simulator

## Introduction

The AUG simulator allows us to simulate the behavior of autonomous underwater gliders so as to analyse their performance and motions, both in a two dimensional as well as three dimensional setting. A PID controller is also implemented for pitch and rudder control.

## Installation

```txt
# clone the project
git clone https://github.com/Bhaswanth-A/AUG-Simulator.git

# install dependencies
pip install -r requirements.txt

# run simulator
python3 main.py
```

## Usage

```txt
usage: main.py [-h] [-i] [-m MODE] [-c CYCLE] [-g GLIDER] [-a ANGLE]
               [-s SPEED] [-pid PID] [-r RUDDER] [-sr SETRUDDER]
               [-p [PLOT ...]]

An Autonomous Underwater Glider Simulator.

options:
  -h, --help            show this help message and exit
  -i, --info            give full information in each cycle
  -m MODE, --mode MODE  set mode as 2D, 3D, or waypoint
  -c CYCLE, --cycle CYCLE
                        number of desired cycles in sawtooth trajectory
  -g GLIDER, --glider GLIDER
                        desired glider model ['slocum']
  -a ANGLE, --angle ANGLE
                        desired glider angle
  -s SPEED, --speed SPEED
                        desired glider speed
  -pid PID, --pid PID   enable or disable PID pitch control
  -r RUDDER, --rudder RUDDER
                        enable or disable rudder
  -sr SETRUDDER, --setrudder SETRUDDER
                        desired rudder angle. Defaults to 10 degrees
  -p [PLOT ...], --plot [PLOT ...]
                        variables to be plotted [3D, all, x, y, z, omega1,
                        omega2, omega3, vel, v1, v2, v3, rp1, rp2, rp3, mb,
                        phi, theta, psi]


```

## TO-Do
- [x] Vertical plane simulations
- [x] 3D simulations
- [x] Params in YAML/py file
- [x] Rudder implementation
- [x] PID pitch control
- [x] PID heading control
- [x] Waypoint following


