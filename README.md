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
# to be updates

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
  -c CYCLE, --cycle CYCLE
                        number of desired cycles in sawtooth trajectory
  -g GLIDER, --glider GLIDER
                        desired glider model
  -a ANGLE, --angle ANGLE
                        desired glider angle
  -s SPEED, --speed SPEED
                        desired glider speed
  -p [PLOT ...], --plot [PLOT ...]
                        variables to be plotted [all, x, z, omega, v, rp1,
                        rp3, mb, theta]
```

## TO-Do
- [x] Vertical plane simulations
- [ ] 3D simulations
- [ ] PID
- [ ] Export sim data
- [x] Params in YAML/py file
- [ ] Incorporate "An underwater glider flight simulator" paper
- [ ] Myring equations?
- [ ] rw


