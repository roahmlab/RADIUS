# RADIUS

## Installing Other Dependecies
* Install Docker Engine, following the instructions on the [Docker website](https://docs.docker.com/desktop/install/linux-install/).
* Download the following into the top-level `util` directory:
  * `coinhsl-2021.05.05.tar.gz` from [here](https://www.hsl.rl.ac.uk/ipopt/) (ACADEMIC LICENCE)
  * `CUDA_Highway_frs.mat`, `CUDA_LeftTurn_frs.mat`, `my_const.mat`,`mycolormap2.mat` from the data folder [here](https://drive.google.com/drive/folders/1ibX50vBhmrv0MuBMZBl0nQztU6D2rbf3?usp=share_link)

## Running the Simulation
* In `Full_Size_Vehicle_Simulation/`, run `./build-docker.sh` to build the Docker
* In `Full_Size_Vehicle_Simulation/simulator/`, run `./run-docker.sh` to run the Docker
* Inside docker in the `/simulator` directory, run `./run-matlab.sh` to invoke MATLAB. MATLAB needs to be activated at the first time when it runs in the docker. In case MATLAB requires account verification but fails to automatically open a browser, try to log in your account and verify the account [online](https://matlab.mathworks.com/) by following steps below.
  1. Click your name initial on the right top corner.
  2. Click 'My Account'.
  3. Click your license number.
  4. Click 'Install and Activate', then click 'View Current Activations' under 'RELATED TASKS'.
  5. Click 'Activate a Computer', then fill in the activation form. Note to choose 'R2022a' as the release version. 
* Open another terminal, and get the 'CONTAINER ID' using `docker ps`. Run `docker exec -it replace_with_your_CONTAINER_ID /bin/bash` to enter the same running docker container, and run `./ninja-cpp-opt.sh` to build the necessary MEX file.
* Run [fully_preprocess_frs.m](https://github.com/roahmlab/RADIUS/blob/main/Full_Size_Vehicle_Simulation/simulator/online/fully_preprocess_frs.m) to generate the C++ condensed version of the FRSes needed. Note that setting the ```is_left_turn``` variable to ```false``` and ```true``` processes FRSes for the highway simulations and left turning simulations respectively.
* In MATLAB invoked by `./run-matlab.sh`, open [RADIUS_simulation_main.m](https://github.com/roahmlab/RADIUS/blob/main/Full_Size_Vehicle_Simulation/simulator/online/RADIUS_simulation_main.m) and run the script. Note that setting the ```scene_type``` variable to 1 and 2 runs the highway simulations and left turning simulations respectively.

