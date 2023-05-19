# RADIUS: Risk-Aware, real-time trajectory Design In Uncertain Scenarios
RADIUS is a risk-aware real-time trajectory planning framework for autonomous driving. 
RADIUS uses Offline zonotope-based reachability analysis on the full order vheicle dynamics to compute the corresponding control-parametrized, over-approximative Forward Reachable Sets (FRS). 
Real-time trajectory planning is achieved by solving an optimization framework in real-time with the pre-computed, control-parametrized FRS being used to ensure vehicle safety up to a given threshold.
The link to the project website can be found [here](https://roahmlab.github.io/Risk_RTD/).

<p align="center">
  <img height="300" src="/util/RADIUS_bigpic.png"/>
</p>

**Authors:** Jinsun Liu* (jinsunl@umich.edu), Challen Enninful Adu* (enninful@umich.edu), Lucas Lymburner (llymburn@umich.edu), Vishrut Kaushik (vishrutk@umich.edu), Lena Trang (ltrang@umich.edu) and Ram Vasudevan (ramv@umich.edu).

*Equal Contribution

All authors are affiliated with the Robotics department and the department of Mechanical Engineering of the University of Michigan, 2505 Hayward Street, Ann Arbor, Michigan, USA.

# Installation Requirements
RADIUS is built on Ubuntu 20.04 with ROS Noetic Distribution, and the algorithms are implemented in MATLAB and C++17. 
<!--RADIUS has the following required dependencies:
- [Docker](https://www.docker.com/) to download simulation package.
- [CORA 2018](https://tumcps.github.io/CORA/) for Forward Reachable Sets representation and computation.-->

# Overview
## 0. Installation
- Clone the RADIUS git repository and run the following from the top level:
```bash
./download-dependencies.sh
```
- Run [install.m](https://github.com/roahmlab/RADIUS/blob/main/install.m).
- In [`split.m`](https://github.com/roahmlab/RADIUS/blob/main/split.m), replace line 20 with ```cd(your_matlab_directory/toolbox/matlab/strfun)``` and line 22 with ```cd('your_CORA2018_directory/global functions/globOptimization')```.
Notice [CORA2018](https://tumcps.github.io/CORA/) is download automatically by ```./download-dependencies.sh```, so it should appear inside [util](https://github.com/roahmlab/RADIUS/tree/main/util).

<!-- - Run dockerfile (need to modify this line once Lucas uploads the file).
- Lunch docker (need to modify this line once Lucas uploads the file).
- Run MATLAB inside docker (need to modify this line once Lucas uploads the file).
- Add required toolboxes in the search path of MATLAB.
- Run [`install.m`](https://github.com/jinsunl/REFINE/blob/main/install.m).
- In [`split.m`](https://github.com/jinsunl/REFINE/blob/main/split.m), replace line 20 with ```cd(your_matlab_directory/toolbox/matlab/strfun)``` and line 22 with ```cd('your_CORA2018_directory/global functions/globOptimization')``` (need to modify this line once Lucas uploads the file). -->
- Continue installation and setup by following instructions in the README for [Full_Size_Vehicle_Simulation](https://github.com/roahmlab/RADIUS/tree/main/Full_Size_Vehicle_Simulation).

## 1. Offline reachability analysis
RADIUS adopts offline reachability analysis from [REFINE](https://github.com/roahmlab/REFINE). 
See [FRS_generation](https://github.com/roahmlab/RADIUS/tree/main/FRS_generation) for details.

## 2. Simulation
RADIUS is evaluated in simulation on a full-size Front-Wheel-Drive vehicle model.
See [Full_Size_Vehicle_Simulation](https://github.com/roahmlab/RADIUS/tree/main/Full_Size_Vehicle_Simulation) for installation and setup details.


<!-- ## 3. Hardware Implementation
REFINE is evaluated in real hardware testing on a All-Wheel-Drive 1/10th race car robot.
See [Rover_Robot_Implementation](https://github.com/jinsunl/REFINE/tree/main/Rover_Robot_Implementation) for details. 
Hardware Demo can be found [here](https://drive.google.com/drive/folders/1FvGHuqIRQpDS5xWRgB30h7exmGTjRyel?usp=sharing). -->



