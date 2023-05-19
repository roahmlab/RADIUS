classdef RADIUS_AgentHelper_Highway < agentHelper
    %% properties
    properties
        % hard reference bounds for parameters
        upper_bd = 10;% y
        lower_bd = 2;% y
        spdlb = 5;% vx
        spdub = 30;%vx
        
        % epsilon for boundary check
        eps = 0.001;
        
        
        HLP;
        S;
        % reference for a single sampling time, needs this so there is a
        % old reference when a new reference comes in
        vx_des = 5; 
        y_des;
        

        plot_flag = 1;
        % reference data for plot
        ref_Z = [];
        proposed_ref_Z = [];
        t_real_start = [];
        t_proposed_start = [];
        FRS_plot_struct;
        
        
        prev_action = -1;
        cur_t0_idx = 1;
        py_idx = 0;

        FRS_helper_handle = struct;
        

        saved_p = [];
        waypt_hist = [];
        p_hist = [];
        FRS_hist = {};
        mirror_hist = [];
        state_hist = [];
        type_manu_hist = [];
        time_hist = [];
        solve_time_hist = [];

        % auxilary variable for debugging
        fixed_thres = [];
        wpt_global = [];
        wpt_takeover = []; % one time waypoint takeover
        v_rel_cric = 10; % [mph]
    end
    %% methods
    methods
        function AH = RADIUS_AgentHelper_Highway(A,FRS_obj,HLP,varargin)
            AH@agentHelper(A,FRS_obj,varargin{:});
            AH.HLP = HLP;
        end

        function [p,tout] = gen_parameter_standalone(AH, world_info, agent_state,waypts)
            % takes a user input parameter and obstacles in world info,
            % find a nearest action that is safe.(may be itself), also
            % return the distance in parameter space and return a replace
            % flag based on if a proper replace action is found

            tic;
            % check if needs to generate a new plan
            if (AH.cur_t0_idx > 1 && AH.prev_action == 2) || (AH.cur_t0_idx > 2 && AH.prev_action == 3)|| AH.prev_action  == 1
                AH.prev_action = -1;
                AH.cur_t0_idx = 1;
            end
            
            % check if a lane change maneuver is still being executed
            if AH.prev_action ~= -1 
                p = [AH.saved_p(1); AH.saved_p(2); AH.cur_t0_idx ;AH.prev_action];
                AH.cur_t0_idx = AH.cur_t0_idx + 1;
                tout = 0;
                return
            end
            
            % get info of world boundary, waypoint, and vehicle state
            world_x_min = world_info.bounds(1);
            world_x_max = world_info.bounds(2);
            world_y_min = world_info.bounds(3);
            world_y_max = world_info.bounds(4);
            world_dyn_obs = get_world_bounds_as_dyn_obs(...
                world_x_min, world_x_max, ...
                world_y_min, world_y_max);
            waypoint_global = waypts(:,1);
            ego_vehicle_state = agent_state(1:6);

            % gather obstacle info
            env_cars = world_info.envCars;
            env_vel = env_cars(:,2)';
            assert(~isempty(AH.A.time))
            env_x0 = env_cars(:,1)' + AH.A.time(end) .* env_vel;
            env_y0 = env_cars(:,3)';
            zer = zeros(size(env_x0));
            env_h0 = zer;
            obs_length = zer + 4.8; % bring in obs footprint
            obs_width = zer + 2.2;
            env_mat = [env_x0;
                       env_y0;
                       env_h0;
                       env_vel;
                       obs_length;
                       obs_width];
            dist_arr = abs(env_x0 - agent_state(1));

            % throw away obstacles if they are too far away from us
            min_dists = min_dist_between_ego_and_dyn_obs(...
                5, ...
                min([agent_state(4)+5, 30.1]), ...
                env_vel, ...
                agent_state(1), ...
                env_x0, ...
                13);
            env_mat(:, or(min_dists > 14, dist_arr > 30*13)) = []; 
            origin_fp = [0; 0; 0; 0; 4.8; 2.2];
            env_mat(:, all(env_mat == origin_fp, 1)) = [];

            % gather static obstacles
            static_obs_idx = env_mat(:, env_mat(4,:) == 0);
            env_mat(:, env_mat(4,:) == 0) = [];
            world_dyn_obs = [world_dyn_obs static_obs_idx];


            dt_sec = 0.01;
            chooseable_dyn_obs_field_name = "dyn_obs";
            chooseable_mu_sigma_field_name = "mu_sigma";


            % generate pdf structures for c++
            sigma_y = 3.7;
            empty_mu_sig_struct = ...
                dyn_obs_to_mu_sigma_struct(0.0, 0.0, 0.0, ...
                0.0, 1.0, 2.2, 0.01, 0.01, 13.0, sigma_y);
            empty_always_risky = [empty_mu_sig_struct];
            empty_always_risky(:) = [];
            chooseable_obs_curr = struct(...
                chooseable_mu_sigma_field_name, [], ...
                chooseable_dyn_obs_field_name, []);
            chooseable_obs_set = [chooseable_obs_curr];
            chooseable_obs_set(1) = [];
            
            % Note our zonotope reacchable sets have various time
            % durations, thus here we generate pdfs over all possible time
            % intervals of various lengths that may accord with one
            % zonotope reachable set. 
            for obs_idx = 1:size(env_mat, 2)
                dyn_obs_curr = env_mat(:, obs_idx);
                vel = dyn_obs_curr(4);
                x0 = dyn_obs_curr(1);
                y0 = dyn_obs_curr(2);
                h0 = dyn_obs_curr(3);
                len = dyn_obs_curr(5);
                width = dyn_obs_curr(6);
                for i = 1:16
                    t_width_sec = dt_sec * i;
                    curr_mu_sigma_struct = dyn_obs_to_mu_sigma_struct(x0, y0, h0, vel, len, width, dt_sec, t_width_sec, 13.0, sigma_y);
                    if i == 1
                        curr_mu_sigmas_chooseable = zeros([0 0]);
                        curr_mu_sigmas_chooseable = [curr_mu_sigma_struct];
                    else
                        curr_mu_sigmas_chooseable(1, i) = curr_mu_sigma_struct;
                    end
                end
                chooseable_obs_curr = struct(chooseable_mu_sigma_field_name, curr_mu_sigmas_chooseable, chooseable_dyn_obs_field_name, dyn_obs_curr);
                if obs_idx == 1
                    chooseable_obs_set = [chooseable_obs_curr];
                else
                    chooseable_obs_set = [chooseable_obs_set chooseable_obs_curr];
                end
            end

            % this is where online planning is solved in c++
            frs_filepath = "/data/cpp_processed_CUDA_FRS_30-Aug-2022_lessFRS_3rnd_grid24_t_fix.DBG";
            risk_cpp_output = RISK_RTD_MEX(ego_vehicle_state, ...
                waypoint_global, ...
                chooseable_obs_set, ...
                world_dyn_obs, ...
                empty_always_risky, ...
                frs_filepath, ...
                AH.fixed_thres);

            % parsing planning result
            if risk_cpp_output(1) == 1 % success
                manu_type = risk_cpp_output(4);
                param_val_unmirrored = risk_cpp_output(2);
                is_mirrored = logical(risk_cpp_output(3));
                u0 = risk_cpp_output(5);
                u0_idx = risk_cpp_output(6) + 2;
                idx0 = risk_cpp_output(7) + 1;
                idx1 = risk_cpp_output(8) + 1;
                if is_mirrored
                    multiplier = -1;
                    param_val = -param_val_unmirrored;
                else
                    multiplier = 1;
                    param_val = param_val_unmirrored;
                end
                if manu_type == 1
                    p = [param_val; 0; 1; manu_type];
                    AH.cur_t0_idx = 1;
                    AH.prev_action = -1;
                else
                    p = [u0; param_val; AH.cur_t0_idx; manu_type];
                    AH.prev_action = manu_type;
                    AH.cur_t0_idx = 2;
                    AH.py_idx = idx1 + (2 * is_mirrored);
                    AH.saved_p = p;
                end
                AH.wpt_global = waypts(1:2);
                
                tout = 0;
                M = AH.zono_full.M_mega{u0_idx};
                type_manu_all = ["Au","dir","lan"];
                type_text = type_manu_all(manu_type);
                FRS = M(type_text); 
                if size(FRS,1) == 1
                    FRS = FRS';
                end

                % info for visualizing FRS
                FRS = FRS{idx0,idx1};
                AH.FRS_plot_struct.p = param_val;
                AH.FRS_plot_struct.type_manu = manu_type;
                AH.FRS_plot_struct.FRS = FRS;
                AH.FRS_plot_struct.agent_state = agent_state;
                AH.FRS_plot_struct.multiplier = multiplier;
                AH.FRS_plot_struct.mirror_flag = is_mirrored;
                AH.FRS_plot_struct.agent_state = agent_state;
                AH.FRS_plot_struct.multiplier = multiplier;
    
                % data to be saved if needed
                AH.FRS_hist{end+1} = FRS;
                AH.mirror_hist = [AH.mirror_hist is_mirrored];
                AH.type_manu_hist = [AH.type_manu_hist manu_type];
                AH.p_hist = [AH.p_hist param_val];
                AH.waypt_hist = [AH.waypt_hist waypoint_global];
                AH.state_hist = [AH.state_hist agent_state];
                AH.time_hist = [AH.time_hist AH.A.time(end)];
            else % not success
                p = [];
                tout = 0;
            end
        end
        function plot_FRS(AH)
            if isstruct(AH.FRS_plot_struct) && ~isempty(AH.FRS_plot_struct.p)
                p = AH.FRS_plot_struct.p;
                type_manu = AH.FRS_plot_struct.type_manu;
                FRS = AH.FRS_plot_struct.FRS;
                mirror_flag = AH.FRS_plot_struct.mirror_flag;
                agent_state = AH.FRS_plot_struct.agent_state;
                multiplier = AH.FRS_plot_struct.multiplier;
                AH.plot_selected_parameter_FRS(p,type_manu,FRS,mirror_flag,agent_state,multiplier);
            end
        end
        function plot_selected_parameter_FRS(AH,p,type_manu,FRS,mirror_flag,agent_state,multiplier)
            if ~isempty(p)
                %clear data and then plot
                AH.FRS_helper_handle.XData = cell(3,1);
                AH.FRS_helper_handle.YData = cell(3,1);
                if type_manu == 1
                    AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[p; 0],[0 1 0],2);
                else
                    AH.plot_zono_collide_sliced(FRS,mirror_flag,agent_state,[agent_state(4);p *multiplier],[0 1 0],2);
                end
            end
        end
        function p = convert_action_to_parameter(AH,action,discrete_flag)
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
                delta_struct.y_des = action(2)*0.59;         %%%%%%  y
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
            p = AH.update_desired_parameters(delta_struct);
            p(2) = delta_struct.y_des; %skipping the previous step y change
            
        end
        function [T, U, Z]=gen_ref(AH, p, ~,agent_state, ~)
            % generate reference trajectory based on parameter and states
            if ~exist('agent_state','var')
                agent_info = AH.get_agent_info();
                agent_state = agent_info.state(:,end);
            end
            u_cur = agent_state(4) ;
            p_u = p(1);
            p_y = p(2);
            t0_idx = p(3);
            t0 = (t0_idx-1)*AH.t_move;
            type_manu = p(4);
            if type_manu == 3
                [T, U,Z] = gaussian_T_parameterized_traj_with_brake(t0,p_y,p_u,u_cur,[],0,1);
            else
                [T, U,Z] = sin_one_hump_parameterized_traj_with_brake(t0,p_y,p_u,u_cur,[],0,1);
            end
        end
        
        %% helper functions
        function parameters=update_desired_parameters(AH,delta_struct)
            parameters = zeros(2,1);
            agent_info = AH.get_agent_info();
            agent_state = agent_info.state(:,end);
            v_cur = agent_state(4) ;
            y_cur = agent_state(2) ;
            
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
            AH.p_hist = [];
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
            for t_idx = 1:length(FRS.vehRS_save) % 1:10/AH.truncating_factor:length(FRS.vehRS_save)
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
                XY = [h.XData(:) h.YData(:)]; % Create Matrix Of Vectors
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
        end

    end
end

