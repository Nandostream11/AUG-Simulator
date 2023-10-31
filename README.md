# Autonomous Underwater Glider Simulator

## Introduction

The AUG simulator allows us to simulate the behavior of autonomous underwater gliders so as to analyse their performance and motions, both in a two dimensional as well as three dimensional setting.

## Requirements

- Python _version_

## Installation

```txt
# clone the project
git clone https://github.com/Bhaswanth-A/AUG-Simulator.git

# install dependencies
# to be updated

# run simulator
python3 main.py
```

## Usage

```txt
usage: main.py [-h] [-i] [-c CYCLE] [-g GLIDER] [-a ANGLE] [-s SPEED]
               [-p [PLOT ...]]

An Autonomous Underwater Glider Simulator.

options:
  -h, --help            show this help message and exit
  -i, --info            give full information in each cycle
  -m MODE, --mode MODE  set mode as 2D or 3D
  -c CYCLE, --cycle CYCLE
                        number of desired cycles in sawtooth trajectory
  -g GLIDER, --glider GLIDER
                        desired glider model ['slocum']
  -a ANGLE, --angle ANGLE
                        desired glider angle
  -s SPEED, --speed SPEED
                        desired glider speed
  -p [PLOT ...], --plot [PLOT ...]
                        variables to be plotted [3D, all, x, y, z, omega1,
                        omega2, omega3, vel, v1, v2, v3, rp1, rp2, rp3, mb, phi,
                        theta, psi]

```

## TO-Do
- [x] Vertical plane simulations
- [x] 3D simulations
- [ ] PID
- [ ] Export sim data
- [x] Params in YAML/py file

