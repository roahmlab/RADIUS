classdef dynamic_car_world < world
    properties
        buffer = 1 ;

        %used for collision checking
        obstacles_unseen
        obstacles_seen
        obstacle_size_bounds = [0.2 0.6];
        obstacle_rotation_bounds = [0,0] ;
        bounds_as_polyline
        
        % plotting
        obstacle_seen_color = [255 255 255]/255 ;
        obstacle_unseen_color = [1 0.6 0.6] ;
        obstacle_history
        obstacle_alpha = 0.5;
        
        obstacle_size = [4.8, 2] ;
        envCars;
        num_cars = 2;
        num_moving_cars = 2;

        o %this is obstacle centered at 0
        o_no_frs
        v_array = [22 24 26 28 30 32];
        obs_FRS
        obstacles_old % for collision check
        
        lanewidth = 4.0
        
        car_safe_dist = 25; % w.goal(1) - car_safe_dist is the maximum location of a static obstale
        car_max_vis_dist =120;
        
        SIM_MAX_DISTANCE = 800 * 1; % maximum initial location of a moving vehicle
        SIM_MIN_CLEAR_DIST = 25; %minimum initial location of any vehicle. 
        

        car_max_spd = 25;
        car_min_spd = 15;
        
        start_line

        % planning horizon times
        t_move_and_failsafe
        t_move

        obstacle_spawn_guard = 1; %flag that ensures that obstacles don't spawn too close to ego vehicle in the beginning
        
        mu_sigma % array to store pdf data
        full_pdfs %set of pdfs for full horizon
        t_pdf %number of seconds we need pdfs for
        Ts %pdf  time step
        num_pdfs %number of pdfs to pass to risk based methods

        %sampled obstacle positions and times to be used for visualization
        dyn_obs_pos
        dyn_obspostime

        simulation_idx

        pre_specify_envCars = 0;
        smooth_obstraj = 1;

        pdf_timestep %pdf  time step
             
    end
    
    methods
        %% constructor
        function W = dynamic_car_world(varargin)

            % set default properties
            start_line = varargin{2}(1);
            bounds = [start_line -start_line 0 12] ;
            W.start_line = -start_line;
            N_obstacles = 0 ;
            
            % parse other args
            W = parse_args(W,'bounds',bounds,'N_obstacles',N_obstacles,varargin{:}) ;
            
            % check for mapping toolbox
            try
                polyxpoly([],[],[],[]) ;
            catch
                error(['Please make sure the Mapping Toolbox is installed',...
                    ' so you can use polyxpoly.'])
            end
        end
        
        %% Initial Setup of simulation environment 
        % reset function just resets time, setup function should reset everything

        function setup(W,scenario_idx)
            if exist('scenario_idx','var')
                % to see progress
                scenario_rng = RandStream('mt19937ar','Seed',scenario_idx);
            end

            %spawn cars
            W.placeCars(scenario_rng);
            W.obstacle_history = [];

            % get world bounds
            B = W.bounds ;
            xlo = B(1) ; xhi = B(2) ; ylo = B(3) ; yhi = B(4) ;
            W.bounds_as_polyline = make_box([xhi-xlo,yhi-ylo]) + repmat(mean([xhi xlo ; yhi ylo],2),1,5) ;
    
            
            % generate start position on left side of world with initial
            % heading of 0, and make sure it's not too close to the walls
            b = W.buffer ;
            
            xlo = xlo + 2*b ;
            xhi = xhi - 2*b ;
            ylo = ylo + 2*b ;
            yhi = yhi - 2*b ;
            
            %check if start point has not been set and if not generate the start point 
            if isempty(W.start)
                s = [xlo ;
                    rand_range(ylo, yhi) ;
                    0 ] ;
                W.start = s ;
            end
            
            % generate goal position on right side of room
            if isempty(W.goal)
                g = [xhi - 10 ;
                    (ylo+yhi)/2] ;
                W.goal = g ;
            end

            % Check W.envCars struct which contains initial obstacle positions and velocities 
            % and use that to generate initial obstacle locations to be used for collision checking.
            if ~isempty(W.envCars)
                W.N_obstacles = size(W.envCars,1)-1;

                % generate obstacles around room
                N_obs = W.N_obstacles ;
                
                if N_obs > 0
                    
                    O = nan(2, 6*N_obs) ; % preallocate obstacle matrix
                 
                    % specify obstacle size
                    l = W.obstacle_size(1);
                    w = W.obstacle_size(2);
                    
                    o = [-l/2  l/2 l/2 -l/2 -l/2 ;
                        -w/2 -w/2 w/2  w/2 -w/2 ] ;
                    W.o_no_frs = o;

                    obs_count = 2; %start with second obstacle, since in W.envCars the first element is always the ego vehicle

                    %loop through the O object and populate the obstacle vertices which are used for collision checking
                    for idx = 1:6:(6*N_obs-1)

                        if idx < 6*(W.num_cars - W.num_moving_cars -1)
                            c = [W.envCars(obs_count,1); W.envCars(obs_count,3)];
                        else
                            c = W.dyn_obs_pos{obs_count - (W.num_cars - W.num_moving_cars)}(1,:)';
                        end

                        O(:,idx:idx+4) = W.o_no_frs + repmat(c,1,5) ;
                        obs_count = obs_count + 1;
                    end
                    
                    W.obstacles = O ;
                    W.obstacles_old = O ;
                    
                    W.N_obstacles = N_obs ;
                end
                
                W.obstacles_unseen = [] ;
                W.obstacles_seen = W.obstacles;

                % set up plot data
                W.plot_data.obstacles_seen = W.obstacles ;
                W.plot_data.obstacles_unseen = [] ;
                W.plot_data.start = [] ;
                W.plot_data.goal = [] ;
                W.plot_data.goal_zone = [] ;
                W.plot_data.bounds = [] ;

                %reset time
                W.reset();
            end
        end
        

        function placeCars(W,scenario_rng)
            %this function automatically generates the initial positions and velocities of all the obstacle vehicles
            
            if ~W.pre_specify_envCars
                W.envCars = zeros(W.num_cars, 6);
                W.envCars(1,5) = 1; % ego vehicle starts from the first lane
            end

            W.t_pdf=13; %number of seconds we need pdfs for

            W.pdf_timestep = 0.01; %timestep size for the use of pdfs
            W.num_pdfs=W.t_pdf/W.pdf_timestep; 
            W.mu_sigma=cell(1,2);

            W.mu_sigma{1}=zeros(W.num_pdfs*W.num_moving_cars,6); %cell array for pdfs spanning 0.01s
            W.mu_sigma{2}=zeros((W.num_pdfs-1)*W.num_moving_cars,6); %cell array for pdfs spanning 0.02s

            %dynamic obstacle positions and corresponding times
            W.dyn_obs_pos = cell(1,W.num_moving_cars);
            W.dyn_obspostime = cell(1,W.num_moving_cars);

            W.full_pdfs=cell(2,W.num_moving_cars); %full stack of predicted pdfs for each obstacle

            goal_xpos_buff = W.goal(1) + 300; %variable to hold the goal position buffered by 300m
            jj=1; %index for non static obstacles to use to build pdfs
            
            %load random seed for being able to reproduce 
            rng(W.simulation_idx);


            for i = 2:W.num_cars
                laneOverlap = true;

                %make sure vehicles don't spawn on top of each other
                while laneOverlap
                    laneOverlap = false;
                    is_static_obs = i <= (W.num_cars - W.num_moving_cars);

                    if ~W.pre_specify_envCars
                        if W.obstacle_spawn_guard 
                            if is_static_obs
                                MAX_DIST =  W.goal(1) - W.car_safe_dist; % this is the furthest possible position for a static obstacle
                            else
                                MAX_DIST = W.SIM_MAX_DISTANCE; % this is the furthest possible initial position for a moving obstacle
                            end
                            xPos = (MAX_DIST - W.SIM_MIN_CLEAR_DIST) * rand(scenario_rng,1,1)  + W.SIM_MIN_CLEAR_DIST;% start at least 30m away
                        else
                            MAX_DIST =  abs(1000) - W.car_safe_dist;
                            xPos = MAX_DIST * rand(scenario_rng,1,1)  + 80;% start at least 50
                        end
                        
                        if is_static_obs
                            laneNum = randi(scenario_rng,[0,1])*2; % do not spawn static obstacles in the middle lane
                        else
                            laneNum = randi(scenario_rng,[0,2]); % {0, 1, 2} for bottom, middle and top lane
                        end
                        
                        for j = 1:i-1
                            % check if vehicles are too close to each other in the same lane
                            if  (abs(xPos - W.envCars(j,1)) < W.car_safe_dist) && (laneNum == W.envCars(j,5))
                                laneOverlap = true;
                                break;
                            end
                            if xPos < 80 && laneNum == 0
                                laneOverlap = true;
                                break;
                            end
                        end
                    end
                    if ~laneOverlap % id doesn't overlap
                        if ~W.pre_specify_envCars

                            W.envCars(i,1) = xPos; %x position
                            if is_static_obs
                                W.envCars(i,2) = 0; %longitudinal velocity
                            else
                                W.envCars(i,2) = W.car_min_spd+rand(scenario_rng,1,1) *(W.car_max_spd-W.car_min_spd); %longitudinal velocity 
                            end
                            W.envCars(i,3) = W.lanewidth * laneNum; %y position
                            W.envCars(i,4) = W.car_max_spd - rand(scenario_rng,1,1) *(W.car_max_spd-W.car_min_spd)/3; % JL: envCars(:, 4:6) is not used
                            W.envCars(i,5) = laneNum; % lane index
                            W.envCars(i,6) = 0;
                        end
                        
                        %the following code builds pdfs for each obstacle vehicle that spans over a 13s horizon
                        
                        if ~is_static_obs
                            
                            %find out how many pdfs are needed from the obstacle initial position to the buffered end
                            %position (goal_xpos_buff = 1300m)
                             
                            end_pdfs = round((goal_xpos_buff - W.envCars(i,1))/(W.envCars(i,2)*W.pdf_timestep));

                            %generate all the centers for the pdfs
                            trange_0p1 = (round(linspace(0,end_pdfs-1,end_pdfs))*W.pdf_timestep)+(W.pdf_timestep/2);
                            trange_0p02 = round(linspace(1,end_pdfs-1,end_pdfs-1),4)*0.01;
                            fullpdfs_x0p1 = W.envCars(i,1) + W.envCars(i,2)*trange_0p1;%centers for the pdfs to use with the interval FRS
                            fullpdfs_x0p02  = W.envCars(i,1) + W.envCars(i,2)*trange_0p02;%centers for the pdfs to use with the 0.02s interval FRS
                            
                            
                            %generate the mu's and sigmas of each obstacle
                            pdf_cntr_i = [W.envCars(i,1);W.envCars(i,3)];
                            overapp_box_0p1 = W.generate_obs_boxes(pdf_cntr_i,W.envCars(i,2),W.pdf_timestep); 
                            overapp_box_0p02 = W.generate_obs_boxes(pdf_cntr_i,W.envCars(i,2),0.02); 
                            rot_angle = 0; 
                            pdf_params_0p1 = obs_pdf('gaussian', 'rot_angle', rot_angle, 'obs_box', overapp_box_0p1); 
                            pdf_params_0p02 = obs_pdf('gaussian', 'rot_angle', rot_angle, 'obs_box', overapp_box_0p02); 
                            
                            %build corresponding y positions
                            fullpdfs_y0p1 = repmat(W.envCars(i,3),1,length(fullpdfs_x0p1));
                            fullpdfs_y0p02 = repmat(W.envCars(i,3),1,length(fullpdfs_x0p02));
                            
                            %extract sigmas 
                            sigma_1=[pdf_params_0p1.sigma(1,1),pdf_params_0p1.sigma(1,2),pdf_params_0p1.sigma(2,2)];
                            sigma_2=[pdf_params_0p02.sigma(1,1),pdf_params_0p02.sigma(1,2),pdf_params_0p02.sigma(2,2)];
                            sigmas_0p1 = repmat(sigma_1,length(fullpdfs_x0p1),1);
                            sigmas_0p02 = repmat(sigma_2,length(fullpdfs_x0p02),1);
                            jjs1 = repmat(jj,length(fullpdfs_x0p1),1);
                            jjs2 = repmat(jj,length(fullpdfs_x0p02),1);

                            W.full_pdfs{1,jj}=[fullpdfs_x0p1',fullpdfs_y0p1',sigmas_0p1,jjs1];
                            W.full_pdfs{2,jj}=[fullpdfs_x0p02',fullpdfs_y0p02',sigmas_0p02,jjs2];

                            %Generate random obstacle trajectories.
                            %%%W.smooth_obstraj = 1 makes them smooth trajectories 
                            %%%W.smooth_obstraj = 0 samples obs positions from PDFs
                            sigma = pdf_params_0p1.sigma;
                            mu_cur = W.full_pdfs{1,jj}(1,1:2)'; 
                            c_cur = mvnrnd(mu_cur,sigma);
                            if W.smooth_obstraj
                                samp_idx = 1:(1/W.pdf_timestep):length(fullpdfs_x0p1);
                            else
                                samp_idx = 1:1:length(fullpdfs_x0p1);
                            end

                            samp_idx = [samp_idx,length(fullpdfs_x0p1)];
                            samp_time = (samp_idx-1)*W.pdf_timestep;

                            d_samptime = diff(samp_time);

                            if any(d_samptime == 0)
                                samp_time(d_samptime == 0) = [];
                                samp_idx(d_samptime == 0) = [];
                            end

                            obstacle_pos = zeros(length(samp_idx),2);
                            obstacle_pos(1,:) = c_cur;

                            for kk=2:(length(samp_idx))
                                speed_invalid=1;                          
                                while speed_invalid
                                    if kk <length(samp_idx)
                                        if W.smooth_obstraj
                                            mu_cur = W.full_pdfs{1,jj}((1/W.pdf_timestep)*(kk-1),1:2)'; %sampling every 1 sec
                                        else
                                            mu_cur = W.full_pdfs{1,jj}((kk-1),1:2)'; 
                                        end
                                        c_cur = mvnrnd(mu_cur,sigma);
                                    else
                                        mu_cur = W.full_pdfs{1,jj}(end,1:2)';   
                                        c_cur = mvnrnd(mu_cur,sigma);
                                    end
                                    if round(norm(c_cur(1)-obstacle_pos(kk-1,1)),4) <= 1.6*round(W.car_max_spd,4)
                                        speed_invalid=0;
                                     end
                                end

                                obstacle_pos(kk,:) = c_cur;
                            end
                            if W.smooth_obstraj
                                t_0p1 = linspace(0,samp_time(end),(samp_time(end)/W.pdf_timestep)+1);
                                x_0p1 = interp1(samp_time,obstacle_pos,t_0p1);
                            else
                                t_0p1 = samp_time;
                                x_0p1 = obstacle_pos;
                            end
                          
                            W.dyn_obs_pos{jj} = x_0p1;
                            W.dyn_obspostime{jj} = t_0p1;
                            jj=jj+1;
                        end

                    end
                end
            end
        end

        function world_info = get_world_info(W,agent_info,~)
            %get info about the current state of the obstacles in the world
            %at the beginning of each planning iteration
            W.vdisp('Getting world info!',3)
            
            if nargin > 1
                zcur = agent_info.state(agent_info.position_indices,end)' ;
                zcur = round(zcur,6) ;
                z = unique(zcur,'rows')' ;
                r = agent_info.sensor_radius ;
            else
                W.vdisp('No agent info provided!',2)
                z = W.start ;
                r = 1 ;
            end
            obs_count = 2;N_obs = size(W.envCars,1)-1;% 1 is ego; 
            
            %generate dynamic obstacle locations and predictions of future
            %locations to be used for collision checking
            O = nan(2,6*N_obs);
            O_future = nan(2,6*N_obs);
            c_idx = 1;
            for idx = 1:6:(6*N_obs-1)
                moving_index=W.num_cars-W.num_moving_cars;
                tstart=agent_info.time(end); %time to use to get the starting pdf to feed into risk planner
                tfut = tstart+W.t_move;

                if obs_count <= moving_index
                    c = [W.envCars(obs_count,1); W.envCars(obs_count,3)];
                    O_future(:,idx:idx+4) = W.o_no_frs + repmat(c,1,5);
                    O(:,idx:idx+4) = W.o_no_frs + repmat(c,1,5);
                else

                c_cur = zeros(W.num_moving_cars,2);
                c_fut = zeros(W.num_moving_cars,2);
                
                num_static =W.num_cars - W.num_moving_cars;
                
                for aa = 1:W.num_moving_cars
                    idx_cur = find(round(W.dyn_obspostime{aa},2)==round(tstart,2));
                    idx_fut = find(round(W.dyn_obspostime{aa},2)==round(tfut,2));
                    if isempty(idx_cur)
                        c_cur(aa,:) = W.dyn_obs_pos{aa}(end,:) + [W.envCars(aa+num_static,2)*(tstart-W.dyn_obspostime{aa}(end)),0];
                        c_fut(aa,:) = W.dyn_obs_pos{aa}(end,:) + [W.envCars(aa+num_static,2)*(tstart-W.dyn_obspostime{aa}(end)+W.t_move),0];
                    elseif isempty(idx_fut)
                        c_cur(aa,:) = W.dyn_obs_pos{aa}(idx_cur,:);
                        c_fut(aa,:) = W.dyn_obs_pos{aa}(idx_cur,:) + [W.envCars(aa+num_static,2)*W.t_move,0];
                    else
                        c_fut(aa,:) = W.dyn_obs_pos{aa}(idx_fut,:);
                        c_cur(aa,:) = W.dyn_obs_pos{aa}(idx_cur,:);
                    end
                
                end

                 O_future(:,idx:idx+4) = W.o_no_frs + repmat(c_fut(c_idx,:)',1,5);
                 O(:,idx:idx+4) = W.o_no_frs + repmat(c_cur(c_idx,:)',1,5);

                 c_idx=c_idx+1;
                end

                obs_count = obs_count + 1;
            end

            O_dynamic = cell(0);
            O_dynamic{1} = O; %pos
            O_dynamic{2} = W.envCars(2:end,2);%vel % skip ego
            W.obstacles_old = W.obstacles_seen;
            W.obstacles_seen = O_future;
            W.obstacles = O_future;
            world_info.obstacles = [];
            world_info.dyn_obstacles = O_dynamic;
            world_info.dyn_obstacles_time = [agent_info.time(end) agent_info.time(end)+W.t_move_and_failsafe];
            world_info.bounds = W.bounds ;
            world_info.start = W.start ;
            world_info.goal = W.goal ;
            world_info.dimension = W.dimension ;
            world_info.obs_FRS = W.obs_FRS;

            world_info.num_cars = W.num_cars;
            world_info.num_moving_cars = W.num_moving_cars;

            idx_0p1 = tstart/W.pdf_timestep+1;
            idx_0p02 = tstart/0.01+1;
            
            %generate PDFs for each obstacle's movement over the full planning horizon
            for i=1:W.num_moving_cars
                traj_length=length(W.full_pdfs{1,i});
                take_length = traj_length - idx_0p1;
                if (idx_0p1 + W.num_pdfs) <= traj_length
                    W.mu_sigma{1}((i-1)*W.num_pdfs+1:i*W.num_pdfs,:)=W.full_pdfs{1,i}(idx_0p1:idx_0p1+W.num_pdfs-1,:);
                    W.mu_sigma{2}((i-1)*(W.num_pdfs-1)+1:i*(W.num_pdfs-1),:)=W.full_pdfs{2,i}(idx_0p02:idx_0p02+W.num_pdfs-2,:);
                elseif take_length<0
                    fillerval = [10000 0 0.1 0 0.1];
                    filler1=[repmat(fillerval,W.num_pdfs,1),repmat(i,W.num_pdfs,1)];
                    filler2=[repmat(fillerval,W.num_pdfs-1,1),repmat(i,W.num_pdfs-1,1)];
                    W.mu_sigma{1}((i-1)*W.num_pdfs+1:i*W.num_pdfs,:)=filler1;
                    W.mu_sigma{2}((i-1)*(W.num_pdfs-1)+1:i*(W.num_pdfs-1),:)=filler2;
                else
                    fillerval = [10000 0 0.1 0 0.1];
                    filler1=[repmat(fillerval,W.num_pdfs-take_length-1,1),repmat(i,W.num_pdfs-take_length-1,1)];
                    filler2=[repmat(fillerval,W.num_pdfs-take_length-1,1),repmat(i,W.num_pdfs-take_length-1,1)];
                    W.mu_sigma{1}((i-1)*W.num_pdfs+1:i*W.num_pdfs,:)=[W.full_pdfs{1,i}(idx_0p1:end,:);filler1];
                    W.mu_sigma{2}((i-1)*(W.num_pdfs-1)+1:i*(W.num_pdfs-1),:)=[W.full_pdfs{2,i}(idx_0p02:end,:);filler2];
                end
            end

            world_info.mu_sigma=W.mu_sigma;

            world_info.time = agent_info.time(end);
            world_info.envCars = W.envCars;
        end
        
        
        %% collision check
        
        function out = collision_check(W,agent_info,check_full_traj_flag)
            
            % by default, don't check the full trajectory
            if nargin < 3
                check_full_traj_flag = false ;
            end
            
            % initialize output (innocent until proven guilty)
            out = 0 ;
            
            % extract agent info
            pos_idx = agent_info.position_indices ;
            h_idx = agent_info.heading_index ;
            fp = agent_info.footprint_vertices ;
            Z = agent_info.state ;
            T = agent_info.time ;
            
            % set up obstacles
            O = [W.obstacles_seen, W.bounds_as_polyline,nan(2,1)] ;
            O_old = [W.obstacles_old, W.bounds_as_polyline,nan(2,1)] ;
            
            % if the agent has moved, we need to check if it crashed;
            % alternatively, we could check the entire trajectory for crashes;
            % finally, if it didn't move, then check the full trajectory
            t_start = W.current_time ;
            if check_full_traj_flag
                t_log = true(size(T)) ;
            else
                t_log = T >= t_start ;
            end
            
            if sum(t_log) > 1
                T = T(t_log) ;
                Z = Z(:,t_log) ;
                
                N_chk = 1; % number of checks to make
                
                t_world = T(1) ; % start time of each check
                
                for chk_idx = 1:N_chk
                    % get the time vector to check
                    t_log_idx = T >= t_world ;
                    T_idx = T(t_log_idx) ;
                    Z_idx = Z(:,t_log_idx) ;
                    
                    % create the time vector for interpolation
                    t_chk = (0:0.02:5 ) + t_world ;
                    t_chk = t_chk(t_chk <= T_idx(end)) ;
                    
                    % create the interpolated trajectory
                    Z_chk = match_trajectories(t_chk,T_idx,Z_idx) ;
                    O_idx = 1;
                    while O_idx + 5 <= size(O_old,2)
                        if abs(O(1,O_idx)-O_old(1,O_idx)) > 100
                            O(:,O_idx:O_idx+5) = []; O_old(:,O_idx:O_idx+5) = [];
                        else
                            O_idx = O_idx + 6;
                        end
                    end
                    O_x = match_trajectories(t_chk,[T_idx(1) T_idx(end)],[O_old(1,:)',O(1,:)'])' ;
                    O_y = match_trajectories(t_chk,[T_idx(1) T_idx(end)],[O_old(2,:)',O(2,:)'])' ;

                    % create a copy of the agent footprint, rotated at each
                    % point of the trajectory
                    X = Z_chk(pos_idx,:) ; % xy positions of trajectory
                    N = size(X,2) ;
                    X = repmat(X(:),1,size(fp,2)) ;
                    
                    % X is rep'd for each pt of agent footprint vertices
                    
                    if ~isempty(h_idx) && (length(agent_info.footprint) > 1)
                        % if there is a heading, and the agent is not a circle,
                        % then rotate each footprint contour
                        H = Z_chk(h_idx,:) ;
                        R = rotmat(double(H)) ;
                        F = R*repmat(fp,N,1) + X ;
                    else
                        % otherwise, just place the footprint contour at each
                        % point of the trajectory
                        F = repmat(fp,N,1) + X ;
                    end
                    
                    Fx = [F(1:2:end,:)' ; nan(1,N)] ;
                    Fy = [F(2:2:end,:)' ; nan(1,N)] ;
                    F = [Fx(:)' ; Fy(:)'] ;
                    ci_bool = ones( 1,size(F,2)/6);
                    for chk_time_idx = 1: size(F,2)/6
                        % check if the resulting contour intersects the obstacles
                        F_idx = (6*(chk_time_idx-1)+1):(6*(chk_time_idx));
                        [ci,cyi] =polyxpoly(F(1,F_idx),F(2,F_idx),O_x(chk_time_idx,:),O_y(chk_time_idx,:)) ;
                        if ~isempty(ci)
                            
                            break
                        else
                            ci_bool(chk_time_idx)   = 0;
                        end
                    end
                    
                    if any(ci_bool)
                        out = 1 ;
                    else
                        out = 0;
                    end
                    
                    % increment the time index
                    t_world = t_world + 5 ;
                end
            end
            
            % update the world time index
            W.current_time = agent_info.time(end) ;
        end

        %% Check if agent is at the goal
        function out = goal_check(W,agent_info)
            out = agent_info.position(1,end)>= W.goal(1);
        end
        %% plotting
        function plot(W,agent_info)
            % set up hold if needed
            hold_check = false ;
            if ~ishold
                hold_check = true ;
                hold on
            end
            obs_count = 2;N_obs = size(W.envCars,1)-1;
            O_seen = nan(2,6*N_obs);

            [X_plot, Y_plot] = W.convert_data_for_plots([nan W.obstacles_seen(1,:)],[nan W.obstacles_seen(2,:)]);

            %generates visualization object for plotting
            if check_if_plot_is_available(W,'obstacles_seen')
                W.plot_data.obstacles_seen.XData = X_plot ;
                W.plot_data.obstacles_seen.YData = Y_plot ;
            else
           % replot obstacles every time.
                seen_data = patch(X_plot, Y_plot, W.obstacle_seen_color);
                W.plot_data.obstacles_seen = seen_data ;
            end


            % plot unsensed obstacles
            O_unseen = W.obstacles_unseen ;
            
            if isempty(O_unseen)
                O_unseen = nan(2,1) ;
            end
            
            if check_if_plot_is_available(W,'obstacles_unseen')
                W.plot_data.obstacles_unseen.XData = O_unseen(1,:) ;
                W.plot_data.obstacles_unseen.YData = O_unseen(2,:) ;
            else
                unseen_data = plot(O_unseen(1,:),O_unseen(2,:),'Color',W.obstacle_unseen_color) ;
                W.plot_data.obstacles_unseen = unseen_data ;
            end
            
            % plot start
            s = W.start ;
            
            % plot goal and goal zone
            g = W.goal ;
            
            % plot bounds
            B = W.bounds_as_polyline ;
            
            if hold_check
                hold off ;
            end
        end
        
        function [X, Y] = convert_data_for_plots(W,x,y)
            idn=find(isnan(x));
            Sz=diff(idn)-1;
            Nmax=max(Sz);
            N=numel(Sz);
            X=zeros(Nmax,N);
            Y=X;
            for i=1:N
                X(1:Sz(i),i)=x(idn(i)+1:idn(i+1)-1);
                Y(1:Sz(i),i)=y(idn(i)+1:idn(i+1)-1);
            end
        end

        function overapp_box = generate_obs_boxes(W,pdf_center,vel,Ts)
           %object to be used to generate pdfs

           l = W.obstacle_size(1);
           w = 0.22*3*2; % 3: 3sigma. 2: plus-minus
                    
           o = [-l/2  l/2 l/2 -l/2 -l/2 ;
            -w/2 -w/2 w/2  w/2 -w/2 ] ;
        
           c_range = (vel*Ts)/2;
        
           obs_cent_k = pdf_center - [c_range;0];
           obs_cent_kp1 = pdf_center + [c_range;0];
        
        
           obs_box_k = o + repmat(obs_cent_k,1,5);
           obs_box_kp1 = o + repmat(obs_cent_kp1,1,5);
        
           Lmin = min(obs_box_k(1,:));
           Lmax = max(obs_box_kp1(1,:));
           Wmin = min(obs_box_k(2,:));
           Wmax = max(obs_box_k(2,:));
        
           overapp_box = [Lmin Lmax Lmax Lmin Lmin;
                        Wmin Wmin Wmax Wmax Wmin];

       end

    end
    
end