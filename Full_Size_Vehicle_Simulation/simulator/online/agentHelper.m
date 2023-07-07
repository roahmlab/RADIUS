classdef agentHelper < handle
    %% properties
    properties
        A
        % FRS object
        zono_full
        
        %time info
        t_move
        t_failsafe_move
        
        %saved TUZ for failsafe action when none can be found
        T
        U
        Z
        
        T_hist
        U_hist
        Z_hist
        T_hist_failsafe
        U_hist_failsafe
        Z_hist_failsafe
        
        planned_path
        % flags to pass into children AH
        flags
        
        stuck_count = 0
        verbose = 0 ;
        plot_data
        name = 'agentHelper' ;
    end
    %% methods
    methods
        function AH = agentHelper(A,FRS_obj,varargin)
            % AH = agentHelper(A,FRS_path,varargin) Construct by passing in
            % a agent and FRS file, file handled by children class 
            if nargin > 0
                AH = parse_args(AH,varargin{:}) ;
                AH.A = A;
                AH.zono_full = FRS_obj;%'zono_full_7.13_1spd.mat');
%                 frs_low = load('FRS_Rover_14-Sep-2021_low_spd.mat');
%                 AH.zono_full.M_mega = horzcat(frs_low.M_mega, AH.zono_full.M_mega);
            end
        end
        %% the only useful function in the parent class
        function [action_replaced, replace_distance, stop_sim, k, wp] = advanced_move(AH,action,world_info)
            % [action_replaced, replace_distance, stop_sim, k] = advanced_move(AH,action,world_info)
            % pass in a action from the main file with action range of [-1,1]
            % scale properly and move the agent accordingly, checks for
            % action safety if safety layer is on
%             if AH.A.state(4,end) < 0.3
%                 stop_sim = 1;
%                 action_replaced =1;replace_distance=0; k =[];
%                 return
%             end
            k = AH.convert_action_to_parameter(action,AH.flags.discrete_flag);
            stop_sim = 0;
            if strcmp(AH.flags.safety_layer, 'Z')
                [k, replace_distance, replaced_flag]= AH.adjust(k,world_info);
                if replaced_flag == 0
                    no_replace_action = 0;
                    action_replaced = 0;
                elseif replaced_flag == 1
                    no_replace_action = 0;
                    action_replaced = 1;
                elseif replaced_flag == 2
                    no_replace_action = 1;
                    action_replaced = 1;
                end
            elseif strcmp(AH.flags.safety_layer, 'N')
                replace_distance = 0;
                action_replaced = 0;
                no_replace_action = 0;
            elseif strcmp(AH.flags.safety_layer,'G')
                agent_info = AH.get_agent_info();
                wp = AH.HLP.get_waypoint(world_info,agent_info.state(:,end));
                [time_vec, delta, w_cmd,exit_status] = gen_control_inputs_gpops(AH, world_info, agent_info.state(:,end), wp);
                 AH.A.move(AH.t_move,time_vec', [delta';w_cmd'], [delta';w_cmd']);
                 if exit_status ~= 0
                    stop_sim = 1;
                 end
            elseif strcmp(AH.flags.safety_layer, 'A')
                agent_info = AH.get_agent_info();
                wp = AH.HLP.get_waypoint(world_info,agent_info.state(:,end));
                wp1 = (wp - agent_info.state(1:3,end)) * 0.333 + agent_info.state(1:3,end);
                wp2 = (wp - agent_info.state(1:3,end)) * 0.667 + agent_info.state(1:3,end);
                k = AH.gen_parameter_standalone(world_info,agent_info.state(:,end),[wp,wp1,wp2]);

                replace_distance = 0;
                action_replaced = 1;
                no_replace_action = 0;
                if isempty(k)
                    no_replace_action = 1;
                end
                
            end
            if strcmp(AH.flags.safety_layer, 'A') && ~no_replace_action  
                AH.stuck_count = 0;
                [AH.T, AH.U, AH.Z] = AH.gen_ref(k);
                if isempty(AH.Z) % if ccpba
                    AH.Z = [];
                else
                    AH.Z(3,:) = AH.Z(3,:) + AH.A.state(3, end);
                end
                
                

                %propagate state forward
                AH.A.move(AH.t_move,AH.T, AH.U, AH.Z);

                %for the single merging scenraio, as a temp fix we extend
                %the trajctory by 1s from a 5s lanechange to a 6s
                %lanechange, and append extra FRSs
                if (isequal(agent_info.type,'CCPBA_agent') || isequal(agent_info.type,'Cantelli_agent')) && isprop(AH,'traj_extend') && AH.traj_extend
                    t_add = 0.01:0.01:1;
                    Z_ref = AH.A.state(:,end) + [AH.A.state(4,end)*t_add;zeros(6,length(t_add))];
                    T_ref = t_add;
                    U_ref = zeros(2,length(Z_ref));
                    % don't call the integrator, just assume the agent
                    % perfectly executes the reference trajectory
                    
                    % append the reference trajectory to the agent's
                    % current trajectory
                    AH.A.commit_move_data(T_ref,Z_ref,T_ref,U_ref,Z_ref) ;

                    %update obstacle vehicle locations for collsion
                    %checking to account for 8s braking time
                    tstart = AH.S.t_now;
                    tfut = tstart + AH.t_move + t_add(end);
                    N_obs = size(AH.S.W.envCars,1)-1;% 1 is ego; 
                    num_static = AH.S.W.num_cars-AH.S.W.num_moving_cars;
                    start_pt  = 6*(num_static-1)+1;
                    aa=1;
                    O_future = AH.S.W.obstacles_seen;

                    for idx = start_pt:6:(6*N_obs-1)
                        idx_cur = find(round(AH.S.W.dyn_obspostime{aa},2)==round(tstart,2));
                        idx_fut = find(round(AH.S.W.dyn_obspostime{aa},2)==round(tfut,2));
%                         h = AH.S.W.envCars(aa,4);
                        h=0;
                        rot_h = [cos(h) -sin(h);
                                     sin(h) cos(h)];
                        if isempty(idx_cur)
%                             c_cur(aa,:) = AH.S.W.dyn_obs_pos{aa}(end,:) + [AH.S.W.envCars(aa+num_static,2)*(tstart-AH.S.W.dyn_obspostime{aa}(end)),0];
                            c_fut(aa,:) = AH.S.W.dyn_obs_pos{aa}(end,:) + (rot_h*[AH.S.W.envCars(aa+num_static,2),0]'*(tstart-AH.S.W.dyn_obspostime{aa}(end)+t_add(end)))';
                        elseif isempty(idx_fut)
%                             c_cur(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_cur,:);
                            c_fut(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_cur,:) + (rot_h*[AH.S.W.envCars(aa+num_static,2),0]'*t_add(end))';
                        else
                            c_fut(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_fut,:);
%                             c_cur(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_cur,:);
                        end
                        
                        O_future(:,idx:idx+4) = rot_h*AH.S.W.o_no_frs + repmat(c_fut(aa,:)',1,5);
%                         O(:,idx:idx+4) = AH.S.W.o_no_frs + repmat(c_cur(aa,:)',1,5);
                        aa = aa+1;

                    end
    
                    AH.S.W.obstacles_seen = O_future;
                end

                idx_move = find(AH.t_move > AH.T,1,'last');
                AH.T_hist = [AH.T_hist AH.T(1:idx_move)+AH.A.time(end) - AH.t_move];
                AH.U_hist = [AH.U_hist AH.U(:,1:idx_move)];
                if ~isempty(AH.Z) % if not ccpba
                    AH.Z_hist = [AH.Z_hist AH.Z(:,1:idx_move)];
                end
%                 xlim([AH.A.state(1,end)-18, AH.A.state(1,end)+10])
            elseif strcmp(AH.flags.safety_layer, 'A')
                fprintf('no replacement action found\n');
                AH.stuck_count = AH.stuck_count + 1;
%                 if AH.T(end)-AH.t_move > AH.t_failsafe_move
                    
                if isequal(agent_info.type,'CCPBA_agent')

                    ZZ = AH.updateView(AH.feas); %generate FRS geometric shapes
                    traj = zeros(length(ZZ),2);
                    traj(1,:) = [agent_info.state(1,end), agent_info.state(2,end)];

                    %take center of the FRS at each time step and use
                    %that as the trajectory to follow
                    for ii = 2:length(ZZ)
                        xpt=(max(ZZ{ii}(1,:)) + min(ZZ{ii}(1,:)))/2;
                        ypt=(max(ZZ{ii}(2,:)) + min(ZZ{ii}(2,:)))/2;
                        traj(ii,:) = [xpt,ypt];  
                    end
                    
                    times = linspace(0,length(traj)-1,length(traj))*0.5;

                    trajlocal = traj-traj(1,:); %convert to local frame positions
                    traj = [times',traj]; %add time stamps to states
                    trajlocal = [times',trajlocal]; %add time stamps to states
                    AH.A.ccpba_traj = trajlocal;

                    %fit the an equation for the evolution of y with
                    %respect to x to use to calculate the heading psi
                    ft = fittype('c1*x^5 + c2*x^4 +c3*x^3 + c4*x^2 +c5*x');
                    f_x = fit(AH.A.ccpba_traj(:,2),AH.A.ccpba_traj(:,3),ft,'StartPoint', [0 0 0 0 0]);
                    c_fx = coeffvalues(f_x);
                    x = trajlocal(:,2);
                    psi_traj = atan(5*c_fx(1)*x.^4 + 4*c_fx(2)*x.^3 + 3*c_fx(3)*x.^2 + 2*c_fx(4)*x + c_fx(5));

                    %calculate the longitudinal velocity 
                    xy_points = traj(:,2:3);
                    diffxy = diff(xy_points);

                    distmat = diffxy*diffxy';
                    dist = zeros(length(diffxy)+1,1);
                    for jj=1:length(diffxy)
                        dist(jj) = distmat(jj,jj);
                        if jj>1 && dist(jj) > dist(jj-1)
                            dist(jj) = 0;
                        end
                    end
                    dist = sqrt(dist);
                    ux = dist./0.5;
                    ux(1) = agent_info.state(4,end);
                    zero_dist = find(dist == 0)+1;
                    zero_dist(end) = [];
                    if ~isempty(zero_dist)
                        traj(zero_dist,2:3) = repmat(traj(zero_dist(1)-1,2:3),length(zero_dist),1);
                    end

                    Z_ref = [traj(:,2:3)';psi_traj';ux';zeros(3,length(traj))];
                    U_ref = zeros(2,length(traj));
                    T_ref = times;
                    % don't call the integrator, just assume the agent
                    % perfectly executes the reference trajectory
                    
                    % get the reference trajectory up to time t_move
                    T_loc = 0:AH.A.integrator_time_discretization:times(end) ;
                    Z_loc = match_trajectories(T_loc,T_ref,Z_ref) ;
                    U_loc = match_trajectories(T_loc,T_ref,U_ref) ;
                    
                    % append the reference trajectory to the agent's
                    % current trajectory
                    AH.A.commit_move_data(T_loc,Z_loc,T_loc,U_loc,Z_loc) ;

                    %update obstacle vehicle locations for collsion
                    %checking to account for 8s braking time
                    tstart = AH.S.t_now;
                    tfut = tstart + times(end);
                    N_obs = size(AH.S.W.envCars,1)-1;% 1 is ego; 
                    num_static = AH.S.W.num_cars-AH.S.W.num_moving_cars;
                    start_pt  = 6*(num_static-1)+1;
                    aa=1;
                    O_future = AH.S.W.obstacles_seen;
                    for idx = start_pt:6:(6*N_obs-1)
                        idx_cur = find(round(AH.S.W.dyn_obspostime{aa},2)==round(tstart,2));
                        idx_fut = find(round(AH.S.W.dyn_obspostime{aa},2)==round(tfut,2));
%                         h = AH.S.W.envCars(aa,4);
                        h=0;
                        rot_h = [cos(h) -sin(h);
                                     sin(h) cos(h)];
                        if isempty(idx_cur)
%                             c_cur(aa,:) = AH.S.W.dyn_obs_pos{aa}(end,:) + [AH.S.W.envCars(aa+num_static,2)*(tstart-AH.S.W.dyn_obspostime{aa}(end)),0];
                            c_fut(aa,:) = AH.S.W.dyn_obs_pos{aa}(end,:) + (rot_h*[AH.S.W.envCars(aa+num_static,2),0]'*(tstart-AH.S.W.dyn_obspostime{aa}(end)+times(end)))';
                        elseif isempty(idx_fut)
%                             c_cur(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_cur,:);
                            c_fut(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_cur,:) + (rot_h*[AH.S.W.envCars(aa+num_static,2),0]'*times(end))';
                        else
                            c_fut(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_fut,:);
%                             c_cur(aa,:) = AH.S.W.dyn_obs_pos{aa}(idx_cur,:);
                        end
                        
                        O_future(:,idx:idx+4) = rot_h*AH.S.W.o_no_frs + repmat(c_fut(aa,:)',1,5);
%                         O(:,idx:idx+4) = AH.S.W.o_no_frs + repmat(c_cur(aa,:)',1,5);
                        aa = aa+1;

                    end
    
                    AH.S.W.obstacles_seen = O_future;

                    stop_sim = 1;
                    fprintf('Executing braking maneuver and terminating\n');
                              
                elseif isequal(agent_info.type,'Cantelli_agent')
                    stop_sim = 1;
                    fprintf('Terminating\n');
                else
                    % JL add to match FRS time horizon and ref_traj horizon
                    if ~isempty(AH.FRS_hist)
                        t_itv = interval(AH.FRS_hist{end}.vehRS_save{end}); tub = supremum(t_itv(20));
                        AH.U(:,AH.T>tub) = [];
                        AH.Z(:,AH.T>tub) = [];
                        AH.T(AH.T>tub) = [];

                        idx = find((AH.T-AH.t_move)>0);
                        idx = idx(1)-1;
                        AH.A.move(AH.T(end)-AH.T(idx),AH.T(idx:end)-AH.T(idx), AH.U(:,idx:end), AH.Z(:,idx:end));
                        AH.T = AH.T(idx:end) - AH.T(idx);
                        AH.U = AH.U(:,idx:end);
                        AH.Z = AH.Z(:,idx:end);
                        
                    
                        idx_move = find(AH.t_move > AH.T,1,'last');
                        AH.T_hist_failsafe = [AH.T_hist_failsafe AH.T(1:idx_move)+AH.A.time(end)-AH.t_move];
                        AH.U_hist_failsafe = [AH.U_hist_failsafe AH.U(:,1:idx_move)];
                        AH.Z_hist_failsafe = [AH.Z_hist_failsafe AH.Z(:,1:idx_move)];
                    end
                    stop_sim = 1;
                    fprintf('Terminating\n');

                end
            else
                action_replaced = 1;
                replace_distance = 0;
            end
            
        end
        %% Helper and unimplemented functions
        
        function agent_info = get_agent_info(AH)
            agent_info = AH.A.get_agent_info();
            agent_info.type = class(AH.A);
        end
        function reset(AH,flags,world_start)
            if nargin < 3
                AH.A.reset();
            else
                AH.A.reset(world_start);
            end
            AH.flags = flags;
            error('agent helper reset not implemented')
        end
        function K = convert_action_to_parameter(AH,action)
            error('convert_action_to_parameter not implemented')
        end
        function [K, replace_distance, replaced_flag] = adjust(AH,K,flags)
            error('adjust not implemented')
        end
        function [T, U, Z] = gen_ref(AH,K)
            error('gen_ref not implemented')
        end
        function [k] = gen_parameter(AH,world_info);
            error('gen_parameter not implemented')
        end
        function  plot(AH)
            error('plot not implemented')
        end
        function plot_adjust(AH)
            error('plot_adjust not implemented')
        end
        function plot_A(AH)
            AH.A.plot()
        end

        function vdisp(AH,s,l)
            % Display a string s if the message's verbose level l is greater
            % than or equal to the planner's verbose level.
            if nargin < 3
                l = 1 ;
            end
            if AH.verbose >= l
                if ischar(s)
                    disp(['    AH: ',s])
                else
                    disp('    AH: String not provided!')
                end
            end
        end
    end
end
