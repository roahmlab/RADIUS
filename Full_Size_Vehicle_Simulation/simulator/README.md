# FL_RTD_Highway_Simulation
# New Updates 5/30/22
This repo now hosts comparison between Gpops and FL ZONO.
## How to make nice-looking visulization of simulation
To make plot for a saved simulation summary sim_summary_5-1_11-32-07.902 (sim_summary_(ExitStatus-EpisodeNumber_TimeStamp)),
### Using MATLAB
1. open plot_from_summary.m
2. Change the file name on line 78 to the file name you need
3. Run the file and stop the code half way to adjust the vertical visulizer to desired size
4. Rerun the file and wait for it to complete, saved video is in raw AVI format, thus very large
5. Convert the saved video to mp4 using ffmpeg format for a manageable size
### Using MATLAB & UNREAL
1. open loading_helper_mysim, change the sim_summary name at the top and run the file
2. open ford_viz.slx, double click Simulation 3D Scene Configuration, and click open unreal editor
3. click play in simulink and then click play in unreal, then there should be a animation running. (How many cars show up in simulation depends on how many you set up in loading_helper_mysim)
4. For more help on how to adjust the camera, use [this](https://www.mathworks.com/help/driving/ref/simulation3dsceneconfiguration.html). Press F to unlock the camera (L to turn on lagging), and then use mouse and scroll wheel to drag around for different angle.
5. When you have a good angle, use to F11 to full-screen the game, and use gaming tool WIN+G for recording screen.
### In Premiere pro
1. drag in two videos into the program
2. drag the unreal video on the timeline to set the aspect ratio, and then drag the smaller video on top of the timeline
3. In effect controls, adjust the position of the smaller video to the left of the frame
4. find two time points in the smaller video that match the large video, right click on the video, and adjust the playback speed so the speed of the two videos match
5. Output to mp4
# What does this repository do:
This repository introduces Feedback Linearization method for zonotope RTD.
Zonotope RTD was first introduced in [this](https://arxiv.org/pdf/1904.05728.pdf) paper. Like RTD, Zonotope RTD method combines a planning reachable set(PRS) and a error reachable set(ERS) for a full Forward Reachable Set (FRS) that overapproximates the future states of a dynamic system. Having both ERS and PRS introduce conservatism, but computing FRS directly causes problems because of the complexity of the dynamics. 

Therefore, in this work, we use Feedback Linearization to partially linearize the system dynamics using control. As a result, we simplify the dynamics enough and a fully parameterizable FRS can be computed directly.

We demonstrate the method on a 6 states car system with a simulator environment that our lab built.(Lab simulator not really compatible with highwayAgentHelper anymore) This repository also can run a baseline RTD method. The online portion of the code is also made to be compatible with the Ford simulator found [here](https://github.com/roahmlab/ford_highway_rtd).

Please refer to the Latex doc for full documentation, this documentation is only for running the code.
# Dependencies
All Forward Reachable Sets(FRS) files can be found at Google Driver folder [here](https://drive.google.com/drive/folders/1hb98cht8HXh1EEsOyd3T1uZyS1D91461?usp=sharing), get access from ramv@umich.edu. To get the RTD method's dependencies, please refer to the documentation of this [RTD repo](https://github.com/shaoyifei96/RTD)

[Simulator](https://github.com/skousik/simulator)

[CORA 2018](https://tumcps.github.io/CORA/)

[Zono_RTD_turtlebot_example](https://github.com/pdholmes/zono_RTD_turtlebot_example)

[RTD_FL_Highway_Simulation](https://github.com/shaoyifei96/RTD_FL_Highway_Simulation)

[old RTD](https://github.com/ramvasudevan/RTD)

[RTD_tutorial](https://github.com/skousik/RTD_tutorial)

# Offline Dynamics Generation and FRS generation
## Vehicle Dynamics
***highway_cruising_6_state_agent***: This class inherits from the RTD_agent_2D convention defined in simulator, look in there for how the agent runs during simulation. 

*dynamics* describes the closed loop dynamics model that is integrated online during simulation. It also switches to a low speed model when the when the high speed model becomes inaccurate during braking. The states of the dynamics system are **x,y,h,u,v,r**, which describe the x position, y position, heading in world frame, longitudinal velocity, lateral velocity, and yaw rate in body frame.

***load_const***: This file defines important parameters for the car and should be loaded every time you need some constant defined there. Watch out for naming conflicts!
## Reference Trajectory
        NOTE: Currently the braking trajectory when using symbolic_flag = 1 is not working. The proper symbolic reference for braking is in the **JL addict** section in **veh_dyn_gen_new**

With the dynamics linearized, we now can set a reference trajectory for longitudinal velocity (ud), lateral velocity (vd), and yaw rate (rd). Intuitively, when we need speed change longitudinally. we can change ud and keep other references as zeros. Conversely, if we want lateral motion, we can define a reference in vd and rd, and keep ud constant. We can change these references together to achieve lateral motion and speed change together, as shown in **example_highway_desired_and_braking_trajectory**

We attach a failsafe maneuver to each plan to make sure that each plan comes to a stop, so that we always have a available plan. 

However, To reduce the number of parameter space bins (defined in latex doc), we choose to never input reference input for both speed change and lateral position change together for a given planning time horizon.  All in all, we define three classes of manuevers, called *lane change*, *direction change* and *speed change*. 

Lane Change(6 seconds): **gaussian_T_parameterized_traj_with_brake**: Set reference vd and rd to be the derivative of a Gaussian function, and setting the peak of the reference vd as parameter *Ay*. It also sets a constant ud as parameter value *u0*.

Direction Change(3 seconds): **gaussian_one_hump_parameterized_traj_with_brake**: Set reference vd and rd to be the a Gaussian function, and setting the peak of the reference vd as parameter *Ay*. It also sets a constant ud as parameter value *u0*.

Speed Change(3 seconds):**gaussian_one_hump_parameterized_traj_with_brake**: This maneuver is embedded in lane change. Set reference ud to be the a linear interpolation between initial longitudinal velocity parameter *u0* and desired longitudinal velocity parameter *Au*.

In order to change plan in the middle of a maneuver, we split the reference trajectories into 0.75s intervals, and offset the reference to start 0.75*t0_idx seconds from the start of a plan.

## Computing FRS

**find_proper_reference_brute**: For the lateral motion trajectories, we need to find the proper Ay for a desired amount for lateral movement for each speed for the closed loop system. In this file we use a binary search and simulation to find the proper amount of Ay and save in *dir_change_Ay_info.mat* and *lane_change_Ay_info.mat*.

Since the reference trajectory does not take v0, r0. into account This file also finds

**veh_dyn_gen_new**: contains how we generate the closed loop dynamics in a CORA understandable format. In addition to the parameters for the reference trajectories (u0, Au, Ay), we also add other parameters that depends on the initial conditions (v0, r0). We call Au, Ay the desired parameters, and u0, v0, r0 the initial condition parameters. We place parameter states after the actual states, assign them a uncertainty interval, and give them zero dynamics. During online they can be sliced(defined in latex doc) for improved accuracy.

 **veh_FRS_***: we generate the FRS with the following steps:
 1. For a t0_idx and initial parameters range and desired parameters range, set proper initial condition zonotope.
 2. Calculate the desired maneuver FRS, the braking FRS with high speed model, and the  braking FRS with low speed model.
 3. Data cleaning and saving 

Since uncertainty in initial condition can cause FRS to grow very quickly, we split the range of all possible initial conditions into intervals, and call each one a bin.
   
**examine_FRS**: Before FRS data can be used online, they need to be combined into a accessible format for the online algorithm, this file puts all FRS into a Map structure for fast lookup.
# Online:
We use a receding horizon planning scheme that takes advantage of the FRS computed offline. 

The algorithm tries to find a plan for **AH.t_plan** = 0.75 seconds. If a plan is found, that plan is executed for **AH.t_move** = 0.75 seconds. If we cannot find a next plan in time, or if none of the plans are safe, we execute the failsafe maneuver that is already a part of the previous plan.

Additionally, we assume there is a high level planner(**highway_HLP**) that proposes a waypoint when we need it.

Reiterate with some more deatils:
## Planning Algorithm Steps:
1. Execute a previous unfinished plan for another 0.75 seconds if the remaining of the plan is safe, until the end of the plan.
2. If previous plan is no longer viable, use the centers of the final zonotope in each bin to determine approximately which bin gets the car closest to the waypoint.
3. Optimize within the bin with safety constraints and determine the values of desired parameters.
4. Go to the next best bin if the current bin turns out to be entirely unsafe.
5. If none of the bins contain a safe plan, execute fail-safe.
## Optimization with Safety Constraints
**highwayAgentHelper**: *AH.gen_parameter* is the main planning algorithm proposed above, and *AH.highway_sampling_opt* does the specific continuous optimization for a bin.

We first slice zonotope with respect to initial condition, then represent obstacles as zonotopes and build safety constraints (see latex doc), and lastly optimize over desired parameters using fmincon.

<!-- # Old Doc

[CORA](https://tumcps.github.io/CORA/) checkout to commit 484c54e0d7990312741fddde5a9c9309d3e8808c
[simulator](https://github.com/skousik/simulator)
[RTD_tutorial](https://github.com/skousik/RTD_tutorial)
[RTD](https://github.com/ramvasudevan/RTD) 
[zono_RTD_turtlebot_example](https://github.com/pdholmes/zono_RTD_turtlebot_example)

Use offline_useful folder for final useable code, use offline folder for testing

Notice: 
    Usually code will run load_const.m to obtain a bunch of constants
    Make sure you are not using the same variable names!

Step1:
    example_highway_desired_and_braking_trajectory.m to look at trajecotries.
    2 types of reference trajecotires 
    1) gaussian_T_parameterized_traj_with_brake  do a lane change
    2) gaussian_one_hump_parameterized_traj_with_brake do a direction change
    use find_proper_reference_brute with the right mode to generate proper Ay for each speed range

Step2:
    veh_dyn_gen_lane_change.m generates dynamics files using the gaussian_T_parameterized_traj_with_brake
    This has a tpk of 6 seconds, meant for high speeds
    veh_dyn_gen____________.m generates dynamics files using the gaussian_one_hump_parameterized_traj_with_brake
    This has a tpk of 3 seconds, meant for all speed range with direction change ability.

Step3:
    veh_FRS_lane_change.m generate 6 seconds FRS for lane change above 5m/s
    veh_FRS_dire_change.m generate 3 seconds FRS for direction change
    veh_FRS_sped_change.m generate 3 seconds FRS for spd change

Step4:
    Under Construction: examine_FRS.m (not commented) package all the FRS into a object to be used online -->