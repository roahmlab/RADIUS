# RADIUS

## Installing Other Dependecies
* Install Docker Engine, following the instructions on the [Docker website](https://docs.docker.com/desktop/install/linux-install/).
* Download the following into the top-level `util` directory:
  * `coinhsl-2021.05.05.tar.gz` from [here](https://www.hsl.rl.ac.uk/ipopt/) (ACADEMIC LICENCE)
  * `CUDA_Highway_frs.mat`, `CUDA_LeftTurn_frs.mat`, `my_const.mat`,`mycolormap2.mat` from the data folder [here](https://drive.google.com/drive/folders/1ibX50vBhmrv0MuBMZBl0nQztU6D2rbf3?usp=share_link)

## Running the Simulation
* In `Full_Size_Vehicle_Simulation/`, run `./build-docker.sh` to build the Docker
* In `Full_Size_Vehicle_Simulation/simulator/`, run `./run-docker.sh` to run the Docker
* In the docker in the `/simulator` directory:
  1. Run `./run-matlab.sh`. MATLAB needs to be activated at the first time when it runs in the docker. In case when MATLAB requires for an account verification but fails to automatically open a browser, try to log in your account and verify the account [online](https://matlab.mathworks.com/).
  2. Run `fully_preprocess_frs.m` to generate the C++ condensed version of the FRSes needed
  3. Open another terminal, and get the 'CONTAINER ID' using `docker ps`. Run `docker exec -it replace_with_your_CONTAINER_ID /bin/bash` to enter the same running docker image, and run `./ninja-cpp-opt.sh` to build the necessary MEX file.
  4. In MATLAB invoked by `./run-matlab.sh`, open [RADIUS_simulation_main.m](https://github.com/roahmlab/RADIUS/blob/main/Full_Size_Vehicle_Simulation/simulator/online/RADIUS_simulation_main.m) and run the script. Note that setting the ```scene_type``` variable to 1 and 2 runs the highway simulations and left turning simulation respectively.
