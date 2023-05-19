%% Main Simulation file for RADIUS. Changing the variable scene_type will allow switching between the 3-Lane highway simulation and the Left turning simulation.
%NOTE THAT IF YOU WANT TO SWITCH SCENE TYPES AFTER FIRST RUNNING ON ONECSCENE TYPE YOU SHOULD USE THE "clear all" 
%COMMAND TO ENABLE MATLAB TO LOAD THE RIGHT FRS FOR THE NEW SCENE

close all, clc

%% Simulation parameters: Make changes here
scene_type = 1; % 1: RADIUS 3-lane highway simulation 2: RADIUS Left Turn simulation

plot_fancy_vehicle = 0; %flag to use nice visualization (Slows down simulation visualization due to number of pbject plotted)
plot_pdf = 0; % flag to visualize PDFs (Slows down simulation visualization due to number of objects plotted)
save_video = 0; %flag to save video of visualization
save_result = 0; %flag to save trial data
pre_specify_envCars = 0; %flag to prespecify scenario obstacle and ego locations (used for rerunning particular trials)
use_plot_template = 1;

%set the smooth_obstraj flag to one to generate smooth obstacle vehicle trajectories used for visualization
%setting it to 0 will sample obstacle locations from the PDFs in order to have randomness in obstacle positions and trajectories. 
% This should be set to 0 to reproduce crash rate statistics reported in paper
smooth_obstraj = 0; 

epsilon = 0.05; %risk threshold (0.05 is 5% risk of collision) 

addpath(genpath('/installs/cpp-opt'))
addpath(genpath('/data'))
addpath(genpath('/simulator'))

if scene_type == 1
    frs_filename = 'CUDA_Highway_frs.mat';
    if ~exist('frs','var')
        disp('Loading Highway Driving FRS')
        frs = load(frs_filename) ;
    else
        disp('FRS already loaded') ;
    end
elseif scene_type == 2
    frs_filename = 'CUDA_LeftTurn_frs.mat';
    if ~exist('frs','var')
        disp('Loading Left Turning FRS')
        frs = load(frs_filename) ;
    else
        disp('FRS already loaded') ;
    end
else
    error('Invalid Scene Type selected. Only supported scenes are 1: 3-Lane Highway and 2: Left turning Crossroads.')
end

%% set up required objects for simulation
lanewidth = 4.0;
goal_radius = 6.3;
world_buffer = 1 ; % this is used to make start and goal locations that are not too close to the robot or boundary
verbose_level = 0;

%% automated from here
%setup world environmental parameters
rng(1);
if scene_type == 1
    bounds = [0, 2000, -(lanewidth / 2) - 1, ((lanewidth/2) * 5) + 1];
    t_move = 3;
    t_failsafe_move = 3;
    num_moving_vehicle_list = randsample([5:1:20 25],1000,true); %number of dynamic vehicles in each simulation trial is randomly sampled
elseif scene_type == 2
    bounds = [-40, 40, -40, 40];
    t_move = 4;
    t_failsafe_move = 4;
    num_static_vehicle_list = randsample(4:6,1000,true); %number of static vehicles in each simulation trial is randomly sampled
end

S.eval = 1; %turn on evaluation so summary will be saved to the episode number

hlp_lookahead = 85;

for simulation_idx = 1:10
    for scenario_idx = 1:1000
        num_ego_vehicles = 1;
        assert(num_ego_vehicles == 1, "Cannot have more than one ego vehicle");

        %if you want to manually specify initial car positions for recreating scenarios
        if pre_specify_envCars
            file = dir(pwd);
            for i = 1:length(file)
                if contains(file(i).name, strcat('-',num2str(scenario_idx), '_' ) )
                    load(file(i).name, 'envCars');
                    break
                end
            end
            W.envCars = envCars;
        end


        if scene_type == 1 
            %% 3-Lane Highway Setup
            num_moving_cars = num_moving_vehicle_list(scenario_idx);
            num_static_cars = 3;
            num_total_cars = num_ego_vehicles + num_moving_cars + num_static_cars;

            W = dynamic_car_world('obstacle_spawn_guard', 1, 'bounds', bounds, ...
                'buffer', world_buffer, 'goal', [1010;lanewidth], ...
                'verbose', verbose_level, 'goal_radius', goal_radius, ...
                'num_cars', num_total_cars, 'num_moving_cars', num_moving_cars, ...
                't_move_and_failsafe', t_move+t_failsafe_move, 'simulation_idx', ...
                simulation_idx,'pre_specify_envCars',pre_specify_envCars, ...
                't_move',t_move, 'lanewidth', lanewidth, 'smooth_obstraj',smooth_obstraj); %1 is ego , 2 means 1 obs

            ego_vehicle = highway_cruising_10_state_agent; 
            ego_vehicle. desired_initial_condition = [10;0; 0; 20;0;0;20;0;0;0];
            ego_vehicle.integrator_type= 'ode45';
            ego_vehicle.plot_trajectory_at_time_flag = false;
            HLP = simple_highway_HLP();
            HLP.lookahead = hlp_lookahead; 
            HLP.lanewidth = lanewidth;
            AH = RADIUS_AgentHelper_Highway(ego_vehicle, ...
                frs, ...
                HLP, ...
                't_move', t_move, ...
                't_failsafe_move', t_failsafe_move, ...
                'eps', 0.001, ...
                'verbose', verbose_level);
            AH.fixed_thres = epsilon; % unit: percentage
    
    
            S = rlsimulator(AH,W, 'safety_layer','A', ...
                'plot_fancy_vehicle', plot_fancy_vehicle, 'plot_pdf', plot_pdf,'save_video', save_video, 'save_result', save_result,'threshold', epsilon);  
        elseif scene_type == 2
            %% Left Turning Setup
            num_moving_cars = 4;
            num_static_cars = num_static_vehicle_list(scenario_idx);
            num_total_cars = num_ego_vehicles + num_moving_cars + num_static_cars;
            W = crossroads_world('bounds', bounds, ...
                'buffer', world_buffer, 'goal', [1010;3.7], ...
                'verbose', verbose_level, 'goal_radius', goal_radius, ...
                'num_cars', num_total_cars, 'num_moving_cars', num_moving_cars, ...
                't_move_and_failsafe', t_move+t_failsafe_move, 'simulation_idx', ...
                simulation_idx,'pre_specify_envCars',pre_specify_envCars, ...
                't_move',t_move,'lanewidth', lanewidth, 'smooth_obstraj',smooth_obstraj); %1 is ego , 2 means 1 obs
            
            wpt_takeover = [-25;8;pi];
            ego_vehicle = highway_cruising_10_state_agent;
            ego_vehicle.Kh = 0;
            ego_vehicle.reset([0,-7.5,pi/2,0.1,0,0,0,0,0,0]');
            ego_vehicle.desired_initial_condition = [0,-7.5,pi/2,0.1,0,0,10,0,0,0]';
            ego_vehicle.integrator_type= 'ode45';
            ego_vehicle.plot_trajectory_at_time_flag = false;
            HLP = simple_highway_HLP;
            HLP.lookahead = hlp_lookahead; % 45 for #9, %% LNOTE: lookahead
            AH = RiskAgentHelper_crossroads(ego_vehicle,frs,HLP,'t_move',t_move,'t_failsafe_move',t_failsafe_move,...
                'eps',0.001,'verbose',verbose_level,'wpt_takeover', wpt_takeover);
            
            AH.fixed_thres = epsilon;
            S = rlsimulator_crossroads(AH,W, 'safety_layer','A', ...
                'plot_fancy_vehicle', plot_fancy_vehicle, 'plot_pdf', plot_pdf,'save_video', save_video, 'save_result', save_result,'threshold', epsilon, ...
                'epscur',scenario_idx); 
        end
        
        AH.S = S;        
        if use_plot_template
            plot_template = load('plot_template.mat');
            S.pdf_template = plot_template.pdf_template;
            S.veh_template = plot_template.veh_template;
        end

        S.eval = 1; %turn on evaluation so summary will be saved to the episode number
        
        %% REAL RESET END
    
            rng(scenario_idx+1);
            IsDone4 = 0;
            S.epscur = scenario_idx;
            S.reset();
            for ii = 1:4000
                S.AH.HLP.lookahead = hlp_lookahead;
                AH.planned_path = [linspace(0,1000);repmat([0;0],1,100)];
                switch scene_type
                    case 1
                        [~,~,IsDone,LoggedSignal]=S.step([rand*2-1;rand*2-1]);
                    case 2
                        [~,~,IsDone,LoggedSignal]=S.step([rand*2-1;rand*2-1], wpt_takeover);
                end
                if IsDone == 1 || IsDone == 3 || IsDone == 4 || IsDone == 5
                    %Note: isDone = 3 : Crash
                    %      isDone = 4 : Safely stopped but did not reach goal 
                    %      isDone = 5 : Reached goal!
                    if S.save_video
                        close(S.videoObj);
                    end
                    break
                elseif IsDone == 0
                    break
                end
            end
            pause(1)

    end
end

done = 'Setup Complete'
