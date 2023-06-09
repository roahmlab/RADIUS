# RADIUS: Offline Reachability Analysis

**Note:** offline reachability analysis can also be achieved without the provided [docker image](https://github.com/roahmlab/Risk_RTD/blob/merging-existing/Full_Size_Vehicle_Simulation/Dockerfile) being launched.

## 1. Highway Scenarios
RADIUS utilizes 3 families of desired trajectories to achieve various driving behaviors, namely: *speed changes*, *direction changes* and *lane changes*. 
Each desired trajectory is parametrized by the control parameter `p=[p_u,p_y]` where `p_u` denotes desired longitudinal speed and `p_y` decides desired lateral displacement. 
RADIUS adopts offline reachability analysis for these 3 trajectory families from [REFINE](https://github.com/roahmlab/REFINE) (see [Offline Reachability Analysis](https://github.com/roahmlab/REFINE/tree/main/Offline_Reachability_Analysis) for details).
However, instead of directly using the zonotope reachable sets generated by [REFINE](https://github.com/roahmlab/REFINE), modification is required for RADIUS to evaluate chance constraints using CUDA during online planning.
To generate the necessary data structure for CUDA computation:
- Download [car_frs.mat](https://drive.google.com/drive/folders/1WZbFFhCyhYQlMJxuV4caIzNoa-Q9VZkW) generated by [REFINE](https://github.com/roahmlab/REFINE).
- Run [Highway_FRS2CUDA.m](https://github.com/roahmlab/RADIUS/blob/main/FRS_generation/Highway_FRS2CUDA.m).
 
The resulting reachable sets with data structures for CUDA are provided as [CUDA_Highway_frs.mat](https://drive.google.com/drive/folders/1ibX50vBhmrv0MuBMZBl0nQztU6D2rbf3?usp=sharing).

## 2. Left Turn Scenarios
RADIUS utilizes a new family of desired trajectories to achieve unguarded left turns (see [Appendix B.A](https://arxiv.org/abs/2302.07933) for details).
Closed-loop dynamics of the vehicle model that tracks a left turning maneuver is generated in [veh_dyn_gen_turning.m](https://github.com/roahmlab/RADIUS/blob/main/FRS_generation/left_turning/veh_dyn_gen_turning.m).
To generate the corresponding zonotope reachable sets (CUDA computation considered), run [compute_FRS_turning.m](https://github.com/roahmlab/RADIUS/blob/main/FRS_generation/left_turning/compute_FRS_turning.m).

The resulting reachable sets for left turning are provided [CUDA_LeftTurn_frs.mat](https://drive.google.com/drive/folders/1ibX50vBhmrv0MuBMZBl0nQztU6D2rbf3?usp=sharing).
