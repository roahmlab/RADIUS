classdef RADIUS_AgentHelper_LeftTurn < agentHelper
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
        vx_des = 5; 
        y_des;


        plot_flag = 1;
        % reference data for plot
        ref_Z = [];
        proposed_ref_Z = [];
        t_real_start = [];
        t_proposed_start = [];
        
        
        prev_action = -1;
        cur_t0_idx = 1;


        
        FRS_helper_handle = struct;
 
        saved_p = [];        
        waypt_hist = [];
        p_hist = [];
        FRS_hist = {};
        mirror_hist = [];
        state_hist = [];
        type_manu_hist = [];
        time_hist = [];
        
        FRS_plot_struct;

        fixed_thres = [];

        wpt_global = [];
        wpt_takeover = []; % one time waypoint takeover

    end
    %% methods
    methods
        function AH = RADIUS_AgentHelper_LeftTurn(A,FRS_obj,HLP,varargin)
            AH@agentHelper(A,FRS_obj,varargin{:});
            AH.HLP = HLP;
        end

        function [p,tout] = gen_parameter_standalone(AH, world_info, agent_state,waypts)
            % takes a user input parameter and obstacles in world info,
            % find a nearest action that is safe.(may be itself), also
            % return the distance in parameter space and return a replace
            % flag based on if a proper replace action is found


            if ~isempty(AH.wpt_takeover)
                waypts = AH.wpt_takeover;
            end

            if (AH.cur_t0_idx > 1 && AH.prev_action == 2) || (AH.cur_t0_idx > 2 && AH.prev_action == 3)|| AH.prev_action  == 1 || AH.prev_action == 4
                AH.prev_action = -1;
                AH.cur_t0_idx = 1;
            end
            

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

	    sigma_x_local_fcn = @(local_length, local_width, lanewidth) (local_length/2/3).^2;
	    sigma_y_local_fcn = @(local_length, local_width, lanewidth) (1.32/2/3).^2;
            empty_mu_sig_struct = ...
                dyn_obs_to_mu_sigma_struct(...
                      0.0, 0.0, 0.0, ... % x0, y0, h0
                      0.0, 1.0, 2.2, ... % vel, len, width
                      0.01,          ... % dt_seconds
                      0.01,          ... % t_width_seconds
                      13.0,          ... % t_total_seconds
                      AH.S.W.lanewidth,  ... % lanewidth
                      0.0, ...
                      sigma_x_local_fcn, ...
                      sigma_y_local_fcn  ...
                  );
            always_risky = [empty_mu_sig_struct];
            always_risky(:) = [];

            
            % gather obstacle info
            env_cars = world_info.envCars;
            env_vel = env_cars(:,2)';
            assert(~isempty(AH.A.time))
            env_x0 = env_cars(:,1)';
            env_y0 = env_cars(:,3)';
            zer = zeros(size(env_x0));
            env_h0 = env_cars(:,4)';
            obs_length = zer + 4.8; % bring in obs footprint
            obs_width = zer + 2.2;
            env_mat = [env_x0;
                       env_y0;
                       env_h0;
                       env_vel;
                       obs_length;
                       obs_width];
            static_obs_idx = env_mat(:, env_mat(4,:) == 0);
            env_mat(:, env_mat(4,:) == 0) = [];
            
            dt_sec = 0.01;
            for obs_idx = 1:size(env_mat, 2)
                dyn_obs_curr = env_mat(:, obs_idx);
                vel = dyn_obs_curr(4);
                x0 = dyn_obs_curr(1);
                y0 = dyn_obs_curr(2);
                h0 = dyn_obs_curr(3);
                len = dyn_obs_curr(5);
                width = dyn_obs_curr(6);
                t_now = AH.A.time(end);
                for i = 1:16
                    t_width_sec = dt_sec * i;
                    curr_mu_sigma_struct = dyn_obs_to_mu_sigma_struct(x0, y0, h0, vel, len, width, dt_sec, t_width_sec, 5.5, AH.S.W.lanewidth, t_now, sigma_x_local_fcn, sigma_y_local_fcn);
                    if i == 1
                        curr_mu_sigmas_chooseable = zeros([0 0]);
                        curr_mu_sigmas_chooseable = [curr_mu_sigma_struct];
                    else
                        curr_mu_sigmas_chooseable(1, i) = curr_mu_sigma_struct;
                    end
                end
                chooseable_obs_curr = struct(chooseable_mu_sigma_field_name, curr_mu_sigmas_chooseable, chooseable_dyn_obs_field_name, dyn_obs_curr);
                if obs_idx == 1
                    always_risky = [curr_mu_sigmas_chooseable];
                else
                    always_risky = [always_risky; curr_mu_sigmas_chooseable];
                end
            end

            

            risk_cpp_output = RISK_RTD_MEX(...
                ego_vehicle_state,  ... % EGO STATE GLOBAL FRAME
                waypoint_global,    ... % WAYPOINT GLOBAL FRAME
                chooseable_obs_set, ... % OBSTACLES VARIABLE REFINE/RISK
                world_dyn_obs,      ... % OBSTACLES ALWAYS REFINE
                always_risky,       ... % OBSTACLES ALWAYS RISK
                frs_filepath,       ... % FRS FILE PATH
                AH.fixed_thres,     ... % RISK THRESHOLD
                false,              ... % CHECK MIRRORS?
                false,              ... % USE WAYPOINT MOD. HEURISTC?
                false,              ... % USE FINAL SELECTION HEURISTIC?
                true                ... % USE LEFT COST FUNCTION?
            );


            AH.cur_t0_idx = 1; % table here has extra cols for mirroring
            AH.wpt_global = waypoint_global;
            type_manu = 4; % left turning
            p = [];
            if risk_cpp_output(1) == 1 % success
                param_val_unmirrored = risk_cpp_output(2);
                is_mirrored = logical(risk_cpp_output(3));
                if is_mirrored
                    multiplier = -1;
                    param_val = -param_val_unmirrored;
                else
                    multiplier = 1;
                    param_val = param_val_unmirrored;
                end
               
                % p = [u_final; p_y; t0_idx;type_manu]; 
                p = [AH.zono_full.LeftTurnFRS(1).u_final; param_val; AH.cur_t0_idx; type_manu];
                AH.prev_action = type_manu;
                AH.cur_t0_idx = 2;
                AH.saved_p = p;
                
                % plot info
                FRS = AH.zono_full.LeftTurnFRS(1);
                AH.FRS_plot_struct.p = param_val;
                AH.FRS_plot_struct.type_manu = type_manu;
                AH.FRS_plot_struct.FRS = FRS;
                AH.FRS_plot_struct.agent_state = agent_state;
                AH.FRS_plot_struct.multiplier = multiplier;
                AH.FRS_plot_struct.mirror_flag = is_mirrored;

                % data to be saved if needed
                AH.FRS_hist{end+1} = FRS;
                AH.mirror_hist = [AH.mirror_hist is_mirrored];
                AH.type_manu_hist = [AH.type_manu_hist type_manu];
                AH.p_hist = [AH.p_hist param_val];
                AH.waypt_hist = [AH.waypt_hist waypoint_global];
                AH.state_hist = [AH.state_hist agent_state];
                AH.time_hist = [AH.time_hist AH.A.time(end)];
            end
            tout = 0;
            return;
        end
        
        function plot_FRS(AH)
            if ~isempty(AH.FRS_plot_struct)
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
            p(2) = delta_struct.y_des; %skipping the previous step y change check!!
            
        end
        function [T, U, Z]=gen_ref(AH, p, ~,agent_state, ~)
            % generate reference based on parameter and states
            if ~exist('agent_state','var')
                agent_info = AH.get_agent_info();
                agent_state = agent_info.state(:,end);
            end
            u_cur = agent_state(4) ;
            p_u = p(1); %K = [Au; Ay; t0_idx; type_manu];
            p_y = p(2);
            t0_idx = p(3);
            t0 = (t0_idx-1)*AH.t_move;
            type_manu = p(4);
            if type_manu == 3
                [T, U,Z] = gaussian_T_parameterized_traj_with_brake(t0,p_y,p_u,u_cur,[],0,1);
            elseif type_manu == 4
                [T, U,Z] = parametrized_turnning_with_brake(p_u,p_y,0);
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

            slice_dim = 12;
            k_slice = K(2);

            
            
            for t_idx = 1:length(FRS.vehRS_save) 
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
        end


        
        
    
    end
end
