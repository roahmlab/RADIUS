classdef RiskAgentHelper_crossroads < agentHelper
    %% properties
    properties
        % hard reference bounds for parameters
        upper_bd = 10;% y
        lower_bd = 2;% y
        spdlb = 5;% vx
        spdub = 30;%vx
        
        % epsilon for boundary check
        eps = 0.001;
        
        num_py_lane;
        num_py_dir;
        num_py_left;
        
        HLP;
        % reference for a single sampling time, needs this so there is a
        % old reference when a new reference comes in
        vx_des = 5; 
        y_des;
        
        % from FRS, to determine initial condition and desired condition
        u0_array
        plot_flag = 1;
        pause_flag = 0;
        draw_subplots = 0;
        %get rid of this eventuallylueluelue
        
        
        % reference data for plot
        ref_Z = [];
        FRS_plotting_param = [];
        proposed_ref_Z = [];
        t_real_start = [];
        t_proposed_start = [];
        

        
        prev_action = -1;
        cur_t0_idx = 1;
        py_idx = 0;

        
        t_plan=1.5;
        t_timeout = inf; %600;
        
        fminconopt = [];
        cont_opt = 1;
        constrain_force = 0;
        
        plot_ax;
        
        S 
        
        FRS_helper_handle = struct;
        
        dynamic_obs_plot; %figure 2 dynamics stuff;
        
        %         tmiddle_label
        truncating_factor;
%         zono_dt ;
        
        FRS_cost_function_idx = 150;
        parent_pointer;
        num_times_called = 0;
        saved_K = [];
        
        waypt_hist = [];
        K_hist = [];
        FRS_hist = {};
        mirror_hist = [];
        state_hist = [];
        type_manu_hist = [];
        time_hist = [];
        solve_time_hist = [];
        
        FRS_plot_struct;

        fixed_thres = [];
        wpt_global = [];

        wpt_takeover = []; % one time waypoint takeover

    end
    %% methods
    methods
        function AH = RiskAgentHelper_crossroads(A,FRS_obj,HLP,varargin)
            AH@agentHelper(A,FRS_obj,varargin{:});
            AH.HLP = HLP;
            AH.num_py_left = length(FRS_obj.LeftTurnFRS);
            
            % this is for fmincon. Eeventually will eliminate this once
            % Lucas has the c++ implementation
            o = optimoptions('fmincon') ;
            o.OptimalityTolerance = 1e-3;
            o.MaxIterations = 100 ;
            o.MaxFunctionEvaluations = 100 ;
            o.SpecifyConstraintGradient = true ;
            o.SpecifyObjectiveGradient = true;
            o.Display = 'off';%'iter';
            o.CheckGradients = false;
            o.UseParallel = false;
            AH.fminconopt = o;


            AH.truncating_factor = 1;
        end

        function [K,tout] = gen_parameter_standalone(AH, world_info, agent_state,waypts)
            % takes a user input parameter and obstacles in world info,
            % find a nearest action that is safe.(may be itself), also
            % return the distance in parameter space and return a replace
            % flag based on if a proper replace action is found
%             x_des = waypts(:,1);


            if ~isempty(AH.wpt_takeover)
                waypts = AH.wpt_takeover;
            end

            if (AH.cur_t0_idx > 1 && AH.prev_action == 2) || (AH.cur_t0_idx > 2 && AH.prev_action == 3)|| AH.prev_action  == 1 || AH.prev_action == 4
                AH.prev_action = -1;
                AH.cur_t0_idx = 1;
            end
            
            % world infomation in global frame
            O_all = world_info.obstacles; % static obstacles in the world, NOT static vehicles, so it is []
            dyn_O.obs   = world_info.dyn_obstacles; % world_info.dyn_obstacles{1}: 2-by-(6*obs_num) to indicate current positions of all vehicles. NOTE each vehicle is described by 4+1 vertices and a nan vector, and static vehicle is included.
                                                    % world_info.dyn_obstacles{2}: column vec, velocity of each vehicle
                                                    % world_info.mu_sigma{1}: [mu_sigma_tint1; mu_sigma_tint1; ...] for 13 seconds with each time interval of duration 0.01 sec
                                                    % world_info.mu_sigma{2}: [mu_sigma_tint1; mu_sigma_tint1; ...] for 13 seconds with each time interval of duration 0.02 sec
            dyn_O.num_dyn_cars = world_info.num_moving_cars;
            dyn_O.num_static_cars = world_info.num_cars - world_info.num_moving_cars - 1;


            
            % find the index correesponding to the initial condition (v_ini, h_ini, delta_ini)
            % paramter range and desired parameter range(y_desired, and vd).
            
            % inefficient to generate bounds obstacle everytime
            bounds = world_info.bounds; % boundary of the world: [xmin xmax ymin ymax]
            xlo = bounds(1) ; xhi = bounds(2) ;
            ylo = bounds(3) ; yhi = bounds(4) ;
            
            Blower = [xlo, xhi, xhi, xlo, xlo ; ylo, ylo, ylo-1, ylo-1, ylo] ;
            Bupper = [xlo, xhi, xhi, xlo, xlo ; yhi, yhi, yhi+1, yhi+1, yhi] ;
            B = [Blower, nan(2,1), Bupper, nan(2,1)] ; % make top and bottom bounds into obstacles
            
            O = [O_all B dyn_O.obs{1}(:,1:6*dyn_O.num_static_cars)] ; % all static obstacles including road boundary 
            for i = 1:length(world_info.wall_constraintbox)
                O = [O, world_info.wall_constraintbox{i}, nan(2,1)]; % can't go off the road in a crossroad            
            end

            tsolvetic = tic;
            %find the action to replace starting from current zonotope
            %index
            %% start planning. we only have FRS for p_y>=0, so we have to mirror things for p_y<=0
%             % DO NOT TRY TO UNDERSTAND tb STUFF
%             load my_const.mat
%             %find all options
%             [~,idxu0] = min(abs(AH.u0_array - agent_state(4)));
%             idxAuu0 = idxu0;  
%             M = AH.zono_full.M_mega{idxu0};
%             MAu = AH.zono_full.M_mega{idxAuu0};
%             tbAu = MAu('Autb');                                          %Ay x y h v v r r
%             tbAu_garbage_idx = find(abs(tbAu(1,:) - 5.75)<0.01);
%             tbdir = M('dirtb'); tbdir = repmat(tbdir,[1,2,1]);
%             tblan = M('lantb'); tblan([2,9,12],:) = tblan([2,9,12],:)/2;
%             tblan = repmat(tblan,[1,2,1]);
%             if size(tbdir,1) > 4
%                 tbdir(logical([1 0 1 1 1 1 1 1 0 1 1 0 1 1]),AH.num_py_dir+1:end,:) = -tbdir(logical([1 0 1 1 1 1 1 1 0 1 1 0 1 1]),AH.num_py_dir+1:end,:);
%             else
%                 tbdir(logical([1 0 1 1 ]),AH.num_py_dir+1:end,:) = -tbdir(logical([1 0 1 1 ]),AH.num_py_dir+1:end,:);
%             end
%             if size(tblan,1) > 4
%                 tblan(logical([1 0 1 1 1 1 1 1 0 1 1 0 1 1]),AH.num_py_lane+1:end,:) = -tblan(logical([1 0 1 1 1 1 1 1 0 1 1 0 1 1]),AH.num_py_lane+1:end,:);
%             else
%                 tblan(logical([1 0 1 1]),AH.num_py_lane+1:end,:) = -tblan(logical([1 0 1 1]),AH.num_py_lane+1:end,:);
%             end
%             all_tb = {tbAu, tbdir, tblan};
%             
%        
%             K= [];type_manu_all = ["Au","dir","lan"];
%             
%             
%             if AH.prev_action ~= -1 
%                     K = [AH.saved_K(1); AH.saved_K(2); AH.cur_t0_idx ;AH.prev_action];
%                     AH.cur_t0_idx = AH.cur_t0_idx + 1;
%             end
%             if AH.prev_action == -1
%                 AH.cur_t0_idx = 1; % table here has extra cols for mirroring
%                 AH.wpt_global = waypts(1:2);
% 
%                 %%%%%%%%%%%%%%% JL: deprecated %%%%%%%%%%%%%%%%%%%%%%%%
% %                 [dAu] = rough_manu_cost_estimate(tbAu, agent_state,waypts,AH.cur_t0_idx);
% %                 [ddir] = rough_manu_cost_estimate(tbdir,agent_state,waypts,AH.cur_t0_idx);
% %                 [dlan] = rough_manu_cost_estimate(tblan,agent_state,waypts,AH.cur_t0_idx);
% %                 dAu = dAu + abs(tbAu(1,:) - agent_state(4));
% %                 dAu(tbAu_garbage_idx) = 1000;
% %                   ddir([2,4]) = 1000;
% %                 if abs(agent_state(3))>0ozh.01
% %                     ddir = ddir - 100;
% %                 end
% %                 if min(abs(agent_state(2)-[0,3.7,3.7*2])) >0.2
% %                     dlan = dlan - 100;
% %                 end
% % 
% %                 cannot_move_lateral_speed = 6;
% %                 if agent_state(4) < cannot_move_lateral_speed
% %                     ddir = 1000;
% %                     dlan = 1000;
% %                 end
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 
%                 weight_bias_on_y = 10000;
%                 weight_bias_on_dir = 1.5;
%                 % create dAu
%                 tmp = M('Au');
%                 dAu = zeros(1,length(tmp));
%                 wpt = waypts(1:3,1);
%                 if abs(wpt(2) - agent_state(2)) > 1 && agent_state(4) > 25
%                     wpt(1) = wpt(1) - 30;
%                 end
%                 for i = 1:length(tmp)
%                     xy = [cos(agent_state(3)), -sin(agent_state(3)); sin(agent_state(3)) cos(agent_state(3)) ] * tmp{i}.vehRS_save{tmp{i}.brake_idx1}.Z(1:2,1) + agent_state(1:2);
%                     dAu(i) = sqrt((xy(1) - wpt(1,1))^2 + weight_bias_on_y*(xy(2) - wpt(2,1))^2);
%                 end
%                 wptall = wpt;
% 
%                 % create ddir
%                 tmp = M('dir');
%                 wpt = waypts(1:2,1);
%                 if abs(agent_state(2) - wpt(2)) > 1 || abs(agent_state(3)) > 0.3
%                     wpt(1) = agent_state(1) + agent_state(4)*3;
%                 end
%                 ddir = zeros(1,2*length(tmp));
%                 for i = 1:length(tmp)
%                     xy = [cos(agent_state(3)), -sin(agent_state(3)); sin(agent_state(3)) cos(agent_state(3)) ] * tmp{i}.vehRS_save{tmp{i}.brake_idx1}.Z(1:2,1) + agent_state(1:2);
%                     ddir(i) = weight_bias_on_dir*sqrt((xy(1) - wpt(1,1))^2 + weight_bias_on_y*(xy(2) - wpt(2,1))^2);
%                     xy = [cos(agent_state(3)), -sin(agent_state(3)); sin(agent_state(3)) cos(agent_state(3)) ] * (tmp{i}.vehRS_save{tmp{i}.brake_idx1}.Z(1:2,1) .* [1; -1]) + agent_state(1:2);
%                     ddir(i + length(tmp)) = weight_bias_on_dir*sqrt((xy(1) - wpt(1,1))^2 + weight_bias_on_y*(xy(2) - wpt(2,1))^2);
%                 end
%                 wptall = [wptall, [wpt; 0]];
% 
%                 % create dlan
%                 tmp = M('lan');
%                 wpt = waypts(1:2,1);
%                 if abs(agent_state(2) - wpt(2)) > 1 || abs(agent_state(3)) > 0.3
%                     wpt(1) = agent_state(1) + agent_state(4)*6;
%                 end
%                 dlan = zeros(1,2*size(tmp,1));
%                 for i = 1:size(tmp,1)
%                     xy = [cos(agent_state(3)), -sin(agent_state(3)); sin(agent_state(3)) cos(agent_state(3)) ] * tmp{i}.vehRS_save{tmp{i}.brake_idx1}.Z(1:2,1) + agent_state(1:2);
%                     dlan(i) = sqrt((xy(1) - wpt(1,1))^2 + weight_bias_on_y*(xy(2) - wpt(2,1))^2);
%                     xy = [cos(agent_state(3)), -sin(agent_state(3)); sin(agent_state(3)) cos(agent_state(3)) ] * (tmp{i}.vehRS_save{tmp{i}.brake_idx1}.Z(1:2,1) .* [1; -1]) + agent_state(1:2);
%                     dlan(i+size(tmp,1)) = sqrt((xy(1) - wpt(1,1))^2 + weight_bias_on_y*(xy(2) - wpt(2,1))^2);
%                 end
%                 wptall = [wptall, [wpt; 0]];
% 
% %                 % if velocity is too small, speed up first!
% %                 if agent_state(4) < 15
% %                     ddir(ddir < 1000) = ddir(ddir < 1000) * 0 + 999;
% %                     dlan(dlan < 1000) = dlan(dlan < 1000) * 0 + 999;
% %                 end
% 
% 
%                 
%                 all_dist = {dAu, ddir, dlan};
%                 while any(all_dist{1}<1000) || any(all_dist{2}<1000) || any(all_dist{3}<1000)
%                     [minAud, Auidx] = min(all_dist{1});
%                     [mindird, diridx] = min(all_dist{2});
%                     [minland, lanidx] = min(all_dist{3});
%                     idx_arr = [Auidx diridx lanidx];
%                     [~,type_manu]=min([minAud mindird minland]);
% 
%                     waypts = wptall(:,type_manu);
%                     waypts = [waypts, waypts/3, waypts/3*2];
%                     x_des = waypts(:,1);
% 
%                     
%                     type_text = type_manu_all(type_manu);
%                     if type_manu == 1
%                         FRS = MAu(type_text);
%                     else
%                         FRS = M(type_text);                     %here is original table
%                     end
%                     % online planning next line
%                     disp("manuver "+type_text+" is being solved!!!!!!!!!!!!!!" );
%                     K = AH.highway_sampling_opt( x_des, FRS, all_tb{type_manu}, idx_arr(type_manu), O,dyn_O, world_info.mu_sigma, agent_state, type_manu, AH.cur_t0_idx, waypts);
%                     if isempty(K)
%                         all_dist{type_manu}(idx_arr(type_manu)) = 1000;
%                     else
%                         if type_manu == 1
%                             %Ku, Ky, t0_idx, manu_type
%                             K = [K; 0;1;type_manu];
%                             AH.prev_action = -1;
%                             AH.cur_t0_idx = 1;
%                         else
%                             K = [agent_state(4);K;AH.cur_t0_idx;type_manu];
%                             AH.prev_action = type_manu;
%                             AH.cur_t0_idx = 2;
%                             AH.py_idx = idx_arr(type_manu);
%                             AH.saved_K = K;
%                         end
%                         AH.wpt_global = waypts(1:2);
%                         break;
%                     end
%                 end
%             end

            %% left turn planning starts here

            frs_filepath = "/data/cpp_processed_CUDA_LeftTurn_frs.frs.bin";

            world_dyn_obs = [];
            for wall_cons_idx = 1:length(world_info.wall_constraintbox)
                wall_cons_i = world_info.wall_constraintbox{wall_cons_idx};
                max_vals = max(wall_cons_i');
                min_vals = min(wall_cons_i');
                x_min = min_vals(1);
                x_max = max_vals(1);
                y_min = min_vals(2);
                y_max = min_vals(2);
                world_dyn_obs = [world_dyn_obs, get_world_bounds_as_dyn_obs(...
                                                    x_min, x_max, ...
                                                    y_min, y_max)];
            end
            

            waypoint_global = waypts(:,1);
            ego_vehicle_state = agent_state(1:6);
            chooseable_dyn_obs_field_name = "dyn_obs";
            chooseable_mu_sigma_field_name = "mu_sigma";

            chooseable_obs_curr = struct(...
                chooseable_mu_sigma_field_name, [], ...
                chooseable_dyn_obs_field_name, []);
            chooseable_obs_set = [chooseable_obs_curr];
            chooseable_obs_set(1) = [];

            sigma_y = (1.32/2/3).^2;
            empty_mu_sig_struct = ...
                dyn_obs_to_mu_sigma_struct(0.0, 0.0, 0.0, ...
                0.0, 1.0, 2.2, 0.01, 0.01, 13.0, sigma_y);
            empty_always_risky = [empty_mu_sig_struct];
            empty_always_risky(:) = [];

            

            risk_cpp_output = RISK_RTD_MEX(...
                ego_vehicle_state, ...  % DONE
                waypoint_global, ...    % DONE
                chooseable_obs_set, ... % TODO, EMPTY
                world_dyn_obs, ...      % DONE
                empty_always_risky, ... % TODO, EMPTY
                frs_filepath, ...       % DONE
                1500000000              ); %

            wpt = waypts(:,1);
            AH.cur_t0_idx = 1; % table here has extra cols for mirroring
            AH.wpt_global = waypts(1:2);
            type_manu = 4; % left turning
            for i = 1:AH.num_py_left
                FRS = AH.zono_full.LeftTurnFRS(i);
                
                % generate mu_sigma corresponding to FRS
                if isfield(FRS,'cuda_FRS')
                    t_range = reshape(FRS.cuda_FRS.t_range,2,[]);
                else
                    t_range = [];
                    for frs_idx = 1:length(FRS.vehRS_save)
                        titv = interval(FRS.vehRS_save{frs_idx});
                        titv = titv(20);
                        t_range = [t_range,[infimum(titv); supremum(titv)]];
                    end
                end
                tnow = AH.A.time(end);
                mu_sigma_temp = [];
                for j = 1:size(t_range,2)
                    mu_sigma_temp = [mu_sigma_temp; AH.S.W.gen_mu_sigma_vary_deltaT(tnow,t_range(:,j))];
                end
                mu_sigma = [];
                for j = 1:AH.S.W.num_moving_cars
                    mu_sigma = [mu_sigma;mu_sigma_temp(:,6*j-5:6*j)];
                end
                if ~isempty(mu_sigma)
                    mu_sigma(:,6) = mu_sigma(:,6) - dyn_O.num_static_cars- 1;
                end
                K = AH.left_turn_opt(wpt, FRS, O,dyn_O, mu_sigma, agent_state);
                if ~isempty(K)
                    K = [AH.zono_full.LeftTurnFRS(i).u_final;K;AH.cur_t0_idx;type_manu];
                    AH.prev_action = type_manu;
                    AH.cur_t0_idx = 2;
                    AH.py_idx = i;
                    AH.saved_K = K;
                    AH.wpt_global = waypts(1:2);
                end
            end
            K
            tout = toc(tsolvetic);
        end

        function [K]= left_turn_opt(AH, x_des, FRS, O, dyn_O, mu_sigma, agent_state) 
            K =[];
            mirror_flag = 0;
            multiplier = 1;
            x_des = world_to_local_mirror(agent_state,x_des, mirror_flag);
            if ~isempty(mu_sigma)
                mu_sigma = world_to_local_mirror_pdf(agent_state,mu_sigma, mirror_flag);
            end
            if isempty(O)
                warning('empty obstacles!')
                return
            end
            
            O = world_to_local_mirror(agent_state,O,mirror_flag ); % static obstacle
%             dyn_O.obs{1} = world_to_local_mirror(agent_state,dyn_O.obs{1},mirror_flag ); % dyn_O.obs{1} is not used at all. 

            
            
            if size(FRS,1) == 1
                FRS = FRS';
            end

            itvl = interval(FRS.vehRS_save{1});
            itvl = itvl(12);
            zono_c = mid(itvl);
            zono_g = rad(itvl);
            [K,costval] = AH.find_optimal_cont(O,agent_state,FRS,mirror_flag,4,zono_c,zono_g,x_des, dyn_O, mu_sigma);
            K = multiplier * K;
            
            if ~isempty(K)
                if AH.S.visualize
                    figure(1);
                end
                AH.FRS_plot_struct.k = K;
                AH.FRS_plot_struct.type_manu = 4;
                AH.FRS_plot_struct.FRS = FRS;
                AH.FRS_plot_struct.mirror_flag = mirror_flag;
                AH.FRS_plot_struct.agent_state = agent_state;
                AH.FRS_plot_struct.multiplier = multiplier;
%                 AH.plot_selected_parameter_FRS(K,type_manu,FRS,mirror_flag,agent_state,multiplier);
                AH.waypt_hist = [AH.waypt_hist, x_des];
                AH.K_hist = [AH.K_hist, K];
                AH.FRS_hist{end+1} = FRS;
                AH.mirror_hist = [AH.mirror_hist mirror_flag];
                AH.type_manu_hist = [AH.type_manu_hist 4];
                AH.state_hist = [AH.state_hist agent_state];
                AH.time_hist = [AH.time_hist AH.A.time(end)];
%                 pause(0.1);
            end
            %             end
        end

        
        function [K]= highway_sampling_opt(AH, x_des, FRS,tbdir,diridx_with_mirror,O, dyn_O, mu_sigma, agent_state, type_manu,t0_idx, wayypts) % not used
            K =[];
            if (diridx_with_mirror > AH.num_py_dir  && type_manu == 2) || (diridx_with_mirror > AH.num_py_lane && type_manu ==  3)
                mirror_flag = 1;
                if (diridx_with_mirror > AH.num_py_dir  && type_manu == 2)
                    diridx = diridx_with_mirror - AH.num_py_dir;
                else
                     diridx = diridx_with_mirror - AH.num_py_lane;
                end
                agent_state(5:6) = -agent_state(5:6);
                %                 agent_state(3) = -agent_state(3);
                multiplier = -1;
                
            else
                mirror_flag = 0;
                diridx = diridx_with_mirror;
                multiplier = 1;
            end
            %enable initial check condition r and v if table is extended
            if size(tbdir,1) > 4
                %                 warning('inital condition check disabled')
                
                
                v_lb = tbdir(5,diridx,t0_idx);
                v_ub = tbdir(6,diridx,t0_idx);
                v_eps = abs(v_ub - v_lb) / 10.0;
                if agent_state(5) < v_lb || agent_state(5) > v_ub
                    warning('v0 out of bound')
                    return
                end
                % if previous returned, the following code won't run
                if agent_state(5) < v_lb
                    agent_state(5) = v_lb + v_eps;
                elseif agent_state(5) > v_ub
                    agent_state(5) = v_ub - v_eps;
                end
                
                r_lb = tbdir(7,diridx,t0_idx);
                r_ub = tbdir(8,diridx,t0_idx);
                r_eps = abs(r_ub - r_lb) / 10.0;
                if agent_state(6) < r_lb || agent_state(6) > r_ub
                    warning('r0 out of bound')
                    return
                end
                if agent_state(6) < r_lb
                    agent_state(6) = r_lb + r_eps;
                elseif agent_state(6) > r_ub
                    agent_state(6) = r_ub - r_eps;
                end
            end
            
            %             x_des
            x_des = world_to_local_mirror(agent_state,x_des, mirror_flag);
            wayypts = world_to_local_mirror(agent_state,wayypts, mirror_flag);
            mu_sigma{1} = world_to_local_mirror_pdf(agent_state,mu_sigma{1}, mirror_flag);
            mu_sigma{2} = world_to_local_mirror_pdf(agent_state,mu_sigma{2}, mirror_flag);
            dummy_state = agent_state;
            if isempty(O)
                warning('empty obstacles!')
                return
            end
            


            
            O(1,min(O(1,:))==O(1,:)) = dummy_state(1);
            O(1,max(O(1,:))==O(1,:)) = dummy_state(1)+300;% change x value of the boundary obs 
            O = world_to_local_mirror(dummy_state,O,mirror_flag );
            dyn_O.obs{1} = world_to_local_mirror(dummy_state,dyn_O.obs{1},mirror_flag );

            
            
            if size(FRS,1) == 1
                FRS = FRS';
            end

            FRS = FRS{diridx,t0_idx};
            if AH.cont_opt                                                                    %center value
%                 if size(tbdir,2) == 1
%                     param_gen = 0.05;
%                 else
%                     param_gen = abs((tbdir(1,2,t0_idx)-tbdir(1,1,t0_idx))/2);
%                 end
                itvl = interval(FRS.vehRS_save{1});
                itvl = itvl(11+(type_manu>1));
                zono_c = mid(itvl);
                zono_g = rad(itvl);
                [K,costval] = AH.find_optimal_cont(O,agent_state,FRS,mirror_flag,type_manu,zono_c,zono_g,x_des, dyn_O, mu_sigma);
                K = multiplier * K;
            end
            if ~isempty(K)
                if AH.S.visualize
                    figure(1);
                end
                AH.FRS_plot_struct.k = K;
                AH.FRS_plot_struct.type_manu = type_manu;
                AH.FRS_plot_struct.FRS = FRS;
                AH.FRS_plot_struct.mirror_flag = mirror_flag;
                AH.FRS_plot_struct.agent_state = agent_state;
                AH.FRS_plot_struct.multiplier = multiplier;
%                 AH.plot_selected_parameter_FRS(K,type_manu,FRS,mirror_flag,agent_state,multiplier);
                AH.waypt_hist = [AH.waypt_hist reshape(wayypts,[9,1])];
                AH.K_hist = [AH.K_hist K];
                AH.FRS_hist{end+1} = FRS;
                AH.mirror_hist = [AH.mirror_hist mirror_flag];
                AH.type_manu_hist = [AH.type_manu_hist type_manu];
                AH.state_hist = [AH.state_hist agent_state];
                AH.time_hist = [AH.time_hist AH.A.time(end)];
%                 pause(0.1);
            end
            %             end
        end
        function plot_FRS(AH)
            if ~isempty(AH.FRS_plot_struct)
                k = AH.FRS_plot_struct.k;
                type_manu = AH.FRS_plot_struct.type_manu;
                FRS = AH.FRS_plot_struct.FRS;
                mirror_flag = AH.FRS_plot_struct.mirror_flag;
                agent_state = AH.FRS_plot_struct.agent_state;
                multiplier = AH.FRS_plot_struct.multiplier;
                AH.plot_selected_parameter_FRS(k,type_manu,FRS,mirror_flag,agent_state,multiplier);
            end
        end
        function plot_selected_parameter_FRS(AH,K,type_manu,FRS,mirror_flag,agent_state,multiplier)
            if ~isempty(K)
                %clear data and then plot
                AH.FRS_helper_handle.XData = cell(3,1);
                AH.FRS_helper_handle.YData = cell(3,1);
                if type_manu == 1
                    %                     AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[K; 0],[0 0 1],0);
                    %                     AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[K; 0],[0 0 1],1);
%                     AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[K; 0],[0 0 1],2);
                    AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[K; 0],[0 1 0],2);
                else
                    %                     AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[agent_state(4);K *multiplier],[0 1 0],0);
                    %                     AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[agent_state(4);K *multiplier],[0 1 0],1);
                    AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[agent_state(4);K *multiplier],[0 1 0],2);
                end
            end
        end
        function [c, ceq, gc, gceq] = eval_zono_highway_cons(AH,K,A_con,b_con,s_con,c_k, g_k,start_tic,timeout, slice_idx, GPU_input)
            c = [];
            ceq = [];
            gc = [];
            gceq = [];
            lambdas = (K - c_k)./g_k; % given a parameter, get coefficients on k_slc_G generators
            for k = 1:length(A_con)
                Zk_tmp = A_con{k}*lambdas - b_con{k}; % A*lambda - b <= 0 means inside unsafe set
                [min_c, min_i]=min(-Zk_tmp(1:s_con{k}(1)));
                c = [c min_c]; % this and below can be combined, but lazy so leave them here
                gc= [gc -A_con{k}(min_i,:)'];
                for m = 1: length(s_con{k})-1
                    [min_c, min_i] = min(-Zk_tmp(s_con{k}(m)+1:s_con{k}(m+1)));
                    c = [c min_c];
                    gc = [gc -A_con{k}(min_i+s_con{k}(m),:)'];
                end
            end
            gc = gc / g_k;
            if toc(start_tic) > timeout
                error('Timed out while evaluating constraint function!')
            end
            %we need max(Ax -b)> 0, since fmincon requires nonlin<=0
            %we specify constraint as min(b-Ax)<=0

            if AH.S.W.num_moving_cars ~= 0
                % risk constraint
                [c_risk, ~, gradc_risk, ~] = AH.zonoGPU_risk_con_multizono(K, GPU_input, slice_idx);
                c = [c, c_risk];
                gc = [gc, gradc_risk];
            end
        end
        function [c,ceq,gc,gceq] = zonoGPU_risk_con_multizono(AH, P, GPU_input, slice_idx)
           % IMPORTANT FUNCTION !!!!!!!
            ceq = []; gceq = []; c = []; gc = [];

            mps2mph = 2.23694;
            if isempty(AH.fixed_thres)
                threshold = @(u) 0.01 * (1.53e-6*(u*mps2mph)^4 - 1.84e-4*(u*mps2mph)^3 + 5.49e-3*(u*mps2mph)^2 + 0.011*(u*mps2mph) + 2.23e-3); 
                % the threshold function takes u in [m/s], and 2.23694 transfers [m/s] to [mph]. 
                % This function mono increases before 31.578 [mph], then mono
                % decreases till 59.574 [mph], and starts increasing again.
                dthreshold = @(u) 0.01 * (1.53e-6*(u*mps2mph)^3*4*2.23694 - 1.84e-4*(u*mps2mph)^2*3*mps2mph + 5.49e-3*(u*mps2mph)*2*mps2mph + 0.011*mps2mph);
            else
                threshold = @(u) AH.fixed_thres;
                dthreshold = @(u) 0;
            end


            risk = 0; 
            drisk = 0;
            for i = 1:length(GPU_input.mu_sigma)
                [constr, dconstr, ~, ~] = prob_integration_frs_flzono(GPU_input.num_zono, GPU_input.x0, GPU_input.y0,...
                                                            GPU_input.grid_dx, GPU_input.grid_dy, ...
                                                            GPU_input.block_inzono_list, GPU_input.H1{i}, GPU_input.H2{i}, GPU_input.H4{i},...
                                                            GPU_input.rot_angle, ...
                                                            GPU_input.mu_sigma{i}, GPU_input.cg_p, GPU_input.g_p_x, GPU_input.grid_size, P);
                if any(isnan(constr))
                    i
                end
                risk = risk + constr;
                drisk = drisk + dconstr;
            end
            
            th = 0;
            u0 = GPU_input.u0;
            tracking_error = abs(0); % change this later
            u0range = [u0 - tracking_error, u0 + tracking_error];
            if slice_idx == 11 && P-tracking_error < u0range(1)
                if u0range(2) < 59.574/ mps2mph
                    th1 = threshold(P-tracking_error);
                    th2 = threshold(u0range(2));
                    if th1 < th2
                        c = risk - th1;
                        gc = drisk - dthreshold(P-tracking_error);
                        th = th1;
                    else
                        c = risk - th2;
                        gc = drisk - 0;
                        th = th2;
                    end
                elseif P-tracking_error < 59.574/ mps2mph
                    c = risk - threshold(59.574/ mps2mph);
                    gc = drisk - 0;
                    th =  threshold(59.574/ mps2mph);
                else
                    c = risk - threshold(P-tracking_error);
                    gc = drisk - dthreshold(P-tracking_error);
                    th =  threshold(P-tracking_error);
                end
            elseif slice_idx == 11 && P+tracking_error > u0range(2)
                if P+tracking_error < 59.574/ mps2mph
                    th1 = threshold(Ptracking_error);
                    th2 = threshold(u0range(1));
                    if th1 < th2
                        c = risk - th1;
                        gc = drisk - dthreshold(P+tracking_error);
                        th = th1;
                    else
                        c = risk - th2;
                        gc = drisk - 0;
                        th = th2;
                    end
                elseif u0range(1) < 59.574/ mps2mph
                    c = risk - threshold(59.574/ mps2mph);
                    gc = drisk - 0;
                    th = threshold(59.574/ mps2mph); 
                else
                    c = risk - threshold(u0range(1));
                    gc = drisk - 0;
                    th = threshold(u0range(1));
                end
            else % constant speed maneuver
                if u0range(2) < 59.574/ mps2mph
                    th1 = threshold(u0range(1));
                    th2 = threshold(u0range(2));
                    [~, idx] = min([th1, th2]);
                    c = risk - threshold(u0range(idx));
                    gc = drisk - 0;
                    th = threshold(u0range(idx));
                elseif u0range(1) < 59.574/ mps2mph
                    c = risk - threshold(59.574/ mps2mph);
                    gc = drisk - 0;
                    th = threshold(59.574/ mps2mph);
                else
                    c = risk - threshold(u0range(1));
                    gc = drisk - 0;
                    th = threshold(u0range(1));
                end
            end

            [risk drisk th length(GPU_input.mu_sigma)]
            
        end
        function [avaliable_action_set, costval] = find_optimal_cont(AH,O,agent_state,FRS,mirror_flag,manu_type, zono_c, zono_g,x_des,dyn_O, mu_sigma)
            %             for i = 1:length(FRS.delta_force)
            %                 FRS.delta_force{i} =  deleteAligned(FRS.delta_force{i});
            %             end
            start_tic = tic;
            if ~isfield(FRS,'cuda_FRS')
                FRS.cuda_FRS = [];
            end
            [A_con,b_con,s_con,GPU_input] = AH.generate_constraints(agent_state,O,FRS, FRS.cuda_FRS, mirror_flag,manu_type,dyn_O, mu_sigma);
            
            timeout_t_pk = AH.t_timeout;
            
            if manu_type == 1
                slice_idx = 11;
            else
                slice_idx = 12;
            end

            % CUDA preprocess (ignore CUDA if there's no dynamical obstacles)
            if dyn_O.num_dyn_cars ~= 0
                for i = 1:length(GPU_input.mu_sigma)
                    [GPU_input.x0, GPU_input.y0, GPU_input.H1{i}, GPU_input.H2{i}, GPU_input.H4{i}, ~] = pre_slice_n_IA(GPU_input.num_zono,...
                                                GPU_input.grid_x0, GPU_input.grid_y0, GPU_input.grid_dx, GPU_input.grid_dy, ...
                                                GPU_input.u0v0r0_slice_beta, ...
                                                GPU_input.g_u0_x, GPU_input.g_u0_y, GPU_input.g_v0_x, GPU_input.g_v0_y,...
                                                GPU_input.g_r0_x, GPU_input.g_r0_y,...
                                                GPU_input.block_inzono_list, GPU_input.rot_angle, GPU_input.mu_sigma{i},...
                                                GPU_input.g_p_x, GPU_input.grid_size);
                    if any(isnan(GPU_input.H1{i})) || any(isnan(GPU_input.H2{i})) || any(isnan(GPU_input.H4{i})) 
                        i
                    end
                end
            end
            
            cons = @(K) AH.eval_zono_highway_cons(K, A_con, b_con, s_con,zono_c,zono_g,start_tic, timeout_t_pk, slice_idx, GPU_input) ;
%             cost = @(K) highway_fmincon_cost_fun(K,agent_state, slice_idx,FRS.vehRS_save{AH.FRS_cost_function_idx},x_des);

            if manu_type == 4
                for myidx = length(FRS.vehRS_save):-1:1
                    myc = center(FRS.vehRS_save{myidx});
                    if myc(4)>=FRS.u_final
                        break
                    end
                end
                cost = @(K) highway_fmincon_cost_fun(K,agent_state, slice_idx,FRS.vehRS_save{myidx},x_des);
            else
                cost = @(K) highway_fmincon_cost_fun(K,agent_state, slice_idx,FRS.vehRS_save{FRS.brake_idx1},x_des);
            end
            %             cost = @(K) highway_fmincon_cost_fun(K, x0, x_des);
            lb =[zono_c - zono_g];
            ub =[zono_c + zono_g];
            if manu_type == 1
                ub = min(ub,30); 
            end
            
%             initial_guess = 0.5*(ub+lb);% - AH.eps;

            if manu_type == 1 || manu_type == 4
                initial_guess = zono_c;
            else
                zono_tmp = FRS.vehRS_save{FRS.brake_idx1};
                idx = find(zono_tmp.Z(slice_idx,2:end))+1;
                slice_g = zono_tmp.Z([1,2,slice_idx],idx);
                c = zono_tmp.Z([1,2,slice_idx],1);
                lambda = (x_des(2) - c(2)) / slice_g(2);
                lambda = max(min(lambda,1),-1);
                initial_guess = zono_c + zono_g * lambda;
            end

            try
                [k,fval,exitflag,~] = fmincon(cost, initial_guess, [], [],...
                    [], [], lb, ub, cons,AH.fminconopt) ;
            catch
                exitflag = -1;%timeout
                warning ('optimization times out')
            end
            %
            if exitflag == 1 || exitflag == 2
                avaliable_action_set = k;
                costval = fval;
                toc(start_tic)
            else
                avaliable_action_set = [];
                costval = inf;
            end
        end

        
        function K = convert_action_to_parameter(AH,action,discrete_flag)
            % called in parent class. differnet for each helper
            delta_struct = struct;
            if discrete_flag == true
                if any(newAct == [0 4 8])%M
                    delta_struct.vx_des = 0;
                elseif any(newAct == [1])%A
                    delta_struct.vx_des = 2;
                elseif any(newAct == [2])%B
                    delta_struct.vx_des = -2;
                elseif any(newAct == [3])%H
                    delta_struct.vx_des = -4;
                else
                end
                
                if any(newAct == [4])
                    delta_struct.y_des = 1;
                elseif any(newAct == [8])
                    delta_struct.y_des = -1;
                else
                    delta_struct.y_des = 0;
                end
            else
                delta_struct.vx_des = (action(1)+1)/2*6-4;   %%%%%%  vx
                delta_struct.y_des = action(2)*0.59;%+1)/2*8+2;    %%%%%%  y
            end
            
            if delta_struct.vx_des > 2
                delta_struct.vx_des = 2;
            elseif delta_struct.vx_des < -4
                delta_struct.vx_des = -4;
            end
            
            if delta_struct.y_des > 1
                delta_struct.y_des = 1;
            elseif delta_struct.y_des < -1
                delta_struct.y_des = -1;
            end
            
            %above is for delta vx and delta y, the following is for actual
            %parameters,      : vx and delta y
            K = AH.update_desired_parameters(delta_struct);
            K(2) = delta_struct.y_des; %skipping the previous step y change check!!
            
        end
        function [T, U, Z]=gen_ref(AH, K, real_reference_flag,agent_state, ref_time)
            % generate reference based on parameter and states
            if ~exist('agent_state','var')
                agent_info = AH.get_agent_info();
                agent_state = agent_info.state(:,end);
            end
            if ~exist('ref_time','var')
                ref_time = AH.A.time(end);
            end
            if ~exist('real_reference_flag','var')
                real_reference_flag = 1;
            end
            u_cur = agent_state(4) ;
            %             h_cur = agent_state(3) ;
            y_cur = agent_state(2) ;
            x_cur = agent_state(1) ;
            Au = K(1);%K = [agent_state(4); K;t0_idx;type_manu];
            Ay = K(2);
            t0_idx = K(3);
            
            t0 = (t0_idx-1)*AH.t_move;
            type_manu = K(4);
            load my_const.mat
            %not symbolic mode time argument will be ignored
            %del_y,Au,u0,t,symbolic_flag, scale_Ay_flag                t0
            if type_manu == 3
                [T, U,Z] = gaussian_T_parameterized_traj_with_brake(t0,Ay,Au,u_cur,[],0,1);
            elseif type_manu == 4
                [T, U,Z] = parametrized_turnning_with_brake(Au,Ay,0);
            else
                [T, U,Z] = sin_one_hump_parameterized_traj_with_brake(t0,Ay,Au,u_cur,[],0,1);
            end
            
%             if real_reference_flag
%                 AH.ref_Z=[AH.ref_Z;x_cur+Z(1,:);y_cur+Z(2,:)];% for plotting
%                 AH.t_real_start = [AH.t_real_start;ref_time];
%             else
%                 AH.proposed_ref_Z=[AH.proposed_ref_Z;x_cur+Z(1,:);y_cur+Z(2,:)];% for plotting
%                 AH.t_proposed_start = [AH.t_proposed_start;ref_time];
%             end
            
            
        end
        
        %% helper functions
        function parameters=update_desired_parameters(AH,delta_struct)
            parameters = zeros(2,1);
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            v_cur = agent_state(4) ;
            %             h_cur = agent_state(3) ;
            y_cur = agent_state(2) ;
            %             x_cur = agent_state(1) ;
            %             replaced_flag = 0;
            
            if abs(delta_struct.vx_des)> AH.eps
                AH.vx_des = round(v_cur,1) + delta_struct.vx_des;
            end
            if  AH.vx_des > AH.spdub
                AH.vx_des= AH.spdub;
            end
            if AH.vx_des < AH.spdlb
                AH.vx_des = AH.spdlb;
            end
            if abs(delta_struct.y_des)> AH.eps
                AH.y_des = y_cur + delta_struct.y_des;    % vy_cur
            end
            if  AH.y_des > AH.upper_bd
                AH.y_des= AH.upper_bd;
            end
            if AH.y_des < AH.lower_bd
                AH.y_des = AH.lower_bd;
            end
            parameters(1) = AH.vx_des;
            parameters(2) = AH.y_des-y_cur;
        end
        function plot(AH)
            hold_check = false ;
            if ~ishold
                hold_check = true ;
                hold on
            end
                if ~isempty(AH.planned_path)
                   plot(AH.planned_path(1,:),AH.planned_path(2,:),'k-','LineWidth',1);
%                             plot(AH.proposed_ref_Z(end-1,:),AH.proposed_ref_Z(end,:),'Color','y','LineWidth',3,'LineStyle','--');
                end
            %                 %                 xlim([AH.ref_Z(1,1)-20,AH.ref_Z(1,1)+30]);
            %             end
            %             if ~isempty(AH.ref_Z)
            %                 plot(AH.ref_Z(end-1,:),AH.ref_Z(end,:),'Color',[0 0 0],'LineStyle','-','LineWidth',3);
            %                 %                 plot(AH.ref_Z(end-1,:),AH.ref_Z(end,:),'g--','LineWidth',2);
            %
            %                 %                 xlim([AH.ref_Z(1,1)-20,AH.ref_Z(1,1)+30]);
            %             end
            text(-250,15,"u="+num2str(AH.A.state(4,end))+"m/s",'Color','red','FontSize',15)
            %             text(-250,-20,"v="+num2str(AH.A.state(5,end))+"m/s",'Color','red','FontSize',13)
            
            if hold_check
                hold off ;
            end
            
        end
        
        
        function reset(AH,flags,eps_seed)
            if ~exist('eps_seed','var')
                AH.A.reset();
            else
                rng(eps_seed)
                AH.A.reset();
            end
            AH.y_des = AH.A.state(2,end);
            AH.vx_des = AH.A.state(4,end);
            AH.flags = flags;
            AH.ref_Z = [];
            AH.proposed_ref_Z = [];
            AH.t_real_start = [];
            AH.t_proposed_start = [];
            AH.K_hist = [];
            AH.waypt_hist = [];
            AH.FRS_hist = {};
            AH.mirror_hist = [];
            AH.state_hist = [];
            AH.type_manu_hist = [];
            AH.time_hist = [];
            if ~isempty(AH.HLP)
                AH.HLP.reset();
                if ~isempty(AH.HLP.plot_data.waypoints)
                    AH.HLP.plot_data.waypoints.XData = [];
                    AH.HLP.plot_data.waypoints.YData = [];
                end
                %                 AH.HLP.plot_data.current_waypoint = [];
            end
        end
        
        
        function plot_zono_collide_sliced(AH,FRS,mirror_flag,agent_state,K,color,slice_level)
            %slice level 0, don't slice
            % 1, slice initial condition
            % 2, slice initial and desired
            if K(2) == 0
                slice_dim = 11;
                k_slice = K(1);
            else
                slice_dim = 12;
                k_slice = K(2);
            end
            if AH.draw_subplots
%                 if mirror_flag
%                     figure(2);subplot(3,1,3);
%                 else
%                     figure(2);subplot(3,1,2);
%                 end
                for t_idx = 1:10/AH.truncating_factor:length(FRS.vehRS_save)
                    zono_one = zonotope_slice(FRS.vehRS_save{t_idx}, [7;8;9;slice_dim], [agent_state(4);agent_state(5);agent_state(6);k_slice]);
                    h = plot(zono_one,[1,2],'Color',color);
                end
                
%                 figure(2); subplot(3,1,1);
            end
            
            for t_idx = 1:5:length(FRS.vehRS_save) % 1:10/AH.truncating_factor:length(FRS.vehRS_save)
                if slice_level == 0
                    zono_one = FRS.vehRS_save{t_idx};
                elseif slice_level == 1
                    zono_one = zonotope_slice(FRS.vehRS_save{t_idx}, [7;8;9], [agent_state(4);agent_state(5);agent_state(6)]);
                elseif slice_level == 2
                    zono_one = zonotope_slice(FRS.vehRS_save{t_idx}, [7;8;9;slice_dim], [agent_state(4);agent_state(5);agent_state(6);k_slice]);
                else
                    error('unknown slice_level in plot selected zono');
                end
                h = plot(zono_one,[1,2],'Color',color);
                if mirror_flag
                    h.YData = - h.YData;
                end
                XY = [h.XData(:) h.YData(:)];                                     % Create Matrix Of Vectors
                theta = agent_state(3);
                R=[cos(theta) -sin(theta); sin(theta) cos(theta)]; %CREATE THE MATRIX
                rotXY=XY*R'; %MULTIPLY VECTORS BY THE ROT MATRIX
                Xqr = reshape(rotXY(:,1), size(h.XData,1), []);
                Yqr = reshape(rotXY(:,2), size(h.YData,1), []);
                %SHIFTING
                h.XData = Xqr+agent_state(1);
                h.YData = Yqr+agent_state(2);
                
                AH.FRS_helper_handle.XData{slice_level+1} = [h.YData nan AH.FRS_helper_handle.XData{slice_level+1}];
                AH.FRS_helper_handle.YData{slice_level+1} = [h.XData nan AH.FRS_helper_handle.YData{slice_level+1}];
            end
            
%             figure(2);
        end


        
        function [A_con,b_con,s_con,GPU_input] = generate_constraints(AH, agent_state,O_pts,FRS, cuda_FRS, mirror_flag, manu_type,dyn_O, mu_sigma)
            if ~exist('mirror_flag','var')
                mirror_flag = 0;
            end
            zono_peak_mid_stop  = FRS.vehRS_save;
            n = length(zono_peak_mid_stop);
            
            
            
            O_pts = O_pts';
            obs_dim = [1; 2]; % note that the obstacle exists in the x-y space (not theta or v)
            if manu_type == 1
                p_dim = 11;
                plotc = [0.8,0.8,1];
            else
                if manu_type == 2
                    plotc = [0.8,1,0.8];
                else
                    plotc = [1,0.8,0.8];
                end
                p_dim = 12;
            end
            buffer_dist = 0; % assume no buffer.
            %
            A_con = {};
            b_con = {};
            s_con = {};
            if mirror_flag
                multipler = -1;
            else
                multipler = 1;
            end
%             if AH.plot_flag
%                 if AH.draw_subplots
%                     if mirror_flag
%                         figure(2);subplot(3,1,3)
%                     else
%                         figure(2);subplot(3,1,2)
%                     end
%                     AH.dynamic_obs_plot = plot([0],[0],'m');
%                 end
%             end
    
            GPU_input = cuda_FRS; % just need to add u0v0r0_slice_beta and mu_sigma
            if AH.S.W.num_moving_cars ~= 0
                GPU_input.u0 = agent_state(4);
                GPU_input.u0v0r0_slice_beta = [(agent_state(4) - cuda_FRS.cg_u0v0r0(1)) / cuda_FRS.cg_u0v0r0(4);
                                               (agent_state(5) - cuda_FRS.cg_u0v0r0(2)) / cuda_FRS.cg_u0v0r0(5);
                                               (agent_state(6) - cuda_FRS.cg_u0v0r0(3)) / cuda_FRS.cg_u0v0r0(6)];
            end
            mps2mph = 2.23694;
            no_risk_veh = [];
            for t_idx = [1:1:n,n] % NOTE: this is only for being fast! Lucas will use all frs
%                 if AH.plot_flag && AH.draw_subplots
%                     if mod(t_idx,ceil(10/AH.truncating_factor)) == 0
%                         zono_one = zonotope_slice(zono_peak_mid_stop{t_idx}, [7;8;9], [agent_state(4);agent_state(5);agent_state(6)]);
%                         plot(zono_one,[1,2],'Color',plotc);
%                         AH.dynamic_obs_plot.XData = 0;
%                         AH.dynamic_obs_plot.YData = 0;
%                     end
%                 end
                zono_one = zono_peak_mid_stop{t_idx};
                
                zono_one = zonotope_slice(zono_one, [7;8;9], [agent_state(4);agent_state(5);agent_state(6)]);
                

                Z = zono_one.Z;
                A_obs_array= []; b_obs_array = [];size_array=[]; size_idx = 0;
                c = Z(obs_dim, 1);
                G = Z(:, 2:end);
                for k_idx = 1:length(p_dim)
                    [~, k_col(k_idx)] = find(G(p_dim(k_idx), :) ~= 0); % find "k-sliceable" generators
                end
                k_slc_G = G(obs_dim, k_col);
                k_no_slc_G = G(obs_dim, :);
                k_no_slc_G(:, k_col) = [];
                agent_heading = multipler * agent_state(3);

                %consider each obstacle as a halfspace
                for obs_idx = 1:(size(O_pts,1)+1)/6
                    one_obs = O_pts((obs_idx-1)*6+1:obs_idx*6-1,:);
                    
                    
                    if max(one_obs(:,1)) < -20 || min(one_obs(:,1)) > 250
                        %                     if(max(one_obs(:,1)) - min(one_obs(:,1)) ) < 100 || max(one_obs(:,2))>0
                        
                        continue;
                    end
                    
%                     obs_zono = local_to_zono_with_h(one_obs,agent_heading);
                    obs_zono = local_to_zono_with_h_modify(one_obs,agent_heading); % JL
                    if AH.plot_flag
                        if AH.draw_subplots
                            if t_idx ==1
                                
                                plot(obs_zono,[1 2],'r');
                            end
                        end
                    end
                    %                     [A_obs, b_obs]=zono_to_Ab(obs_zono,Z);
                    %Slice, Do it now! 6 , 7 can be sliced later, taking values vx_des_soft;y_des_soft
                    %Need to get rid of the extra zonos that have negative velocity
                    obstacle = obs_zono.Z;
                    buff_obstacle_c = [obstacle(:, 1) - c];
                    buff_obstacle_G = [obstacle(:, 2:end), k_no_slc_G, buffer_dist*eye(2)]; % obstacle is "buffered" by non-k-sliceable part of FRS
                    buff_obstacle_G(:, ~any(buff_obstacle_G)) = []; % delete zero columns of G
                    buff_obstacle = [buff_obstacle_c, buff_obstacle_G];
                    [A_obs, b_obs] = polytope_PH(buff_obstacle); % turn zonotope into polytope
                    A_obs_array = [A_obs_array;A_obs];
                    b_obs_array = [b_obs_array;b_obs];
                    size_idx = size_idx + length(b_obs);
                    size_array = [size_array;size_idx];
                end
                
                
                %% Dynamics Obstacles
                zonoint = interval(zono_one);
                tint = zonoint(length(zonoint));
                tlb = infimum(tint);
                tub = supremum(tint);

                if manu_type ~= 4
                    if abs(tub - tlb - 0.01) > 1e-4
                        which_mu_sigma = 2;
                        time_idx = round((zono_one.Z(20,1)-0.01)/0.01)+1; % see CORA_2018/contSet/@zonotope/enclose.m for why doing this 
                    else
                        which_mu_sigma = 1;
                        time_idx = round(tlb/0.01)+1;
                    end
                else
                    time_idx = t_idx;
                end
                

                for traffic_idx = 1:dyn_O.num_dyn_cars% each car. NEED TO CHANGE THIS BASED ON CHALLEN'S CHANGE IN DYNAMIC CAR WORLD
                    v_rel = abs(agent_state(4) - dyn_O.obs{2}(traffic_idx+dyn_O.num_static_cars)) * mps2mph;

                    if manu_type == 4
                        mu_sigma_to_use = mu_sigma;
                        v_rel_thres = 10000;
                    else
                        mu_sigma_to_use = mu_sigma{which_mu_sigma};
                        v_rel_thres = 10;
                    end
                    mu_sigma_to_use = mu_sigma_to_use(mu_sigma_to_use(:,6)==traffic_idx, 1:5);
                    mu_sigma_to_use = mu_sigma_to_use(time_idx,:);
                    mu = mu_sigma_to_use(1:2)';
                    sigma = reshape(mu_sigma_to_use([3,4,4,5]),2,2);

                    dist_x = abs(mu(1));
                    dist_y = abs(mu(2));


                    if (v_rel >= v_rel_thres) || (dist_x <= agent_state(4)*1 && dist_y <= 0.75) % flzono type constraint
                        no_risk_veh = unique([no_risk_veh,traffic_idx]);

                        if mu(1) <= -200 || mu(1) >= 200
                            continue
                        end

                        % bound mu_sigma_to_use
                        [evec, eval] = eig(sigma);
                        evec = det(evec)*evec;
                        h_sigma = atan2(evec(2,1),evec(1,1));
                        bla = diag(eval);
                        obs_box = [cos(h_sigma),-sin(h_sigma); sin(h_sigma), cos(h_sigma)]*make_box([bla(1), bla(2)],[0,0]) + mu;
%                         obs_zono_new = local_to_zono_with_h(obs_box',0);
                        obs_zono_new = local_to_zono_with_h_modify(obs_box',0); % JL
                        bla = [cos(agent_heading), sin(agent_heading); -sin(agent_heading), cos(agent_heading)] * diag([4.8/2+0.01*which_mu_sigma*25/2, 2.2/2+0.01*which_mu_sigma*0.5/2]);
                        obs_zono_new = zonotope([obs_zono_new.Z bla]);

%                         if AH.plot_flag
%                             if AH.draw_subplots
%                                 if mod(t_idx,ceil(10/AH.truncating_factor)) == 0
%                                     %                             plot(obs_zono_new,[1 2],'r');
%                                     poly_vert = polygon(obs_zono_new);
%                                     %                                 if min(poly_vert(1,:))>-100
%                                     AH.dynamic_obs_plot.XData = [poly_vert(1,:) nan(1,1) AH.dynamic_obs_plot.XData] ;
%                                     AH.dynamic_obs_plot.YData = [poly_vert(2,:) nan(1,1) AH.dynamic_obs_plot.YData] ;
%                                 end
%                             end
%                         end

                        obstacle = obs_zono_new.Z;
                        buff_obstacle_c = [obstacle(:, 1) - c];
                        buff_obstacle_G = [obstacle(:, 2:end), k_no_slc_G, buffer_dist*eye(2)]; % obstacle is "buffered" by non-k-sliceable part of FRS
                        buff_obstacle_G(:, ~any(buff_obstacle_G)) = []; % delete zero columns of G
                        buff_obstacle = [buff_obstacle_c, buff_obstacle_G];
                        [A_obs, b_obs] = polytope_PH(buff_obstacle); % turn zonotope into polytope
                        A_obs_array = [A_obs_array;A_obs];
                        b_obs_array = [b_obs_array;b_obs];
                        size_idx = size_idx + length(b_obs);
                        size_array = [size_array;size_idx];
                    end
                end
%                 end
%                 if AH.draw_subplots
%                     if mod(t_idx,ceil(10/AH.truncating_factor)) == 0
%                         drawnow
%                     end
%                 end
                A_con{end+1} = A_obs_array*k_slc_G; % now polytope is over coefficients of k_slc_G generators
                b_con{end+1} = b_obs_array;
                s_con{end+1} = size_array;

            end
            
            GPU_input.mu_sigma = cell(1,dyn_O.num_dyn_cars-length(no_risk_veh));
            risk_veh = setdiff(1:dyn_O.num_dyn_cars, no_risk_veh);
            for t_idx = 1:n
                zono_one = zono_peak_mid_stop{t_idx};
                zonoint = interval(zono_one);
                tint = zonoint(length(zonoint));
                tlb = infimum(tint);
                tub = supremum(tint);
                if manu_type ~= 4
                    if abs(tub - tlb - 0.01) > 1e-4
                        which_mu_sigma = 2;
                        time_idx = round((zono_one.Z(20,1)-0.01)/0.01)+1; % see CORA_2018/contSet/@zonotope/enclose.m for why doing this 
                    else
                        which_mu_sigma = 1;
                        time_idx = round(tlb/0.01)+1;
                    end
                else
                    time_idx = t_idx;
                end
                for traffic_idx = 1:length(risk_veh)
                    if manu_type == 4
                        mu_sigma_to_use = mu_sigma;
                    else
                        mu_sigma_to_use = mu_sigma{which_mu_sigma};
                    end
                    mu_sigma_to_use = mu_sigma_to_use(mu_sigma_to_use(:,6)==risk_veh(traffic_idx), 1:5);
                    mu_sigma_to_use = mu_sigma_to_use(time_idx,:);

                    % JL add to test
%                     temp = mu_sigma_to_use(3)*mu_sigma_to_use(5)-mu_sigma_to_use(4)^2;
%                     mu_sigma_to_use(3:5) = mu_sigma_to_use(3:5) / temp / 10;
%                     mu_sigma_to_use(5) = 0.05;

                    %%%%%%%%%%%%% NOTICE: for single iteration planning only %%%%%%%%%
%                     if t_idx >= FRS.brake_idx1
%                         mu_sigma_to_use(1) = 100000;
%                     end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    

                    GPU_input.mu_sigma{traffic_idx} = [GPU_input.mu_sigma{traffic_idx}, mu_sigma_to_use];
                end
            end

        end
    
    end
end
