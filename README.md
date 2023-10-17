# Autonomous Underwater Glider Simulator

## Introduction

The AUG simulator allows us to simulate the behavior of autonomous underwater gliders so as to analyse their performance and motions, both in a two dimensional as well as three dimensional setting.

## Requirements

- Python _version_

## Installation

```
# clone the project
git clone https://github.com/Bhaswanth-A/AUG-Simulator.git

# install dependencies
# to be updates

# run simulator
python3 main.py
```

## Usage

```md
usage: main.py [-h] [-i] [-c CYCLE] [-g GLIDER] [-a ANGLE] [-s SPEED]
               [-p [PLOT ...]]

An Autonomous Underwater Glider Simulator.

options:
  -h, --help            show this help message and exit
  -i, --info            Give full information in each cycle.
  -c CYCLE, --cycle CYCLE
                        Number of desired cycles in sawtooth trajectory.
  -g GLIDER, --glider GLIDER
                        Desired Glider Model.
  -a ANGLE, --angle ANGLE
                        Desired Glider Angle.
  -s SPEED, --speed SPEED
                        Desired Glider Speed.
  -p [PLOT ...], --plot [PLOT ...]
                        Variable to be plotted. Values accepted are [all, x,
                        z, omega, v, rp1, rp3, mb, theta].
```

## TO-Do

- <input type="checkbox" checked> Vertical plane simulations
- <input type="checkbox"> 3D simulations
- <input type="checkbox"> PID
- <input type="checkbox"> Export sim data
- <input type="checkbox" checked> Params in YAML/py file
- <input type="checkbox"> Incorporate "An underwater glider flight simulator" paper
- <input type="checkbox"> Myring equations?
- <input type="checkbox"> rw


