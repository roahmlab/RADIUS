classdef rlsimulator < handle
    properties
        safety_layer = 'A'; % automatic waypoint 
        discrete_flag = 0;
        replace_action = 0; % 0 for original rl lib, 1 for replace with new action, 2 for with punishing reward b/w replaced and original
        plot_sim_flag = 1; %plot simulation figure 1
        plot_AH_flag = 0;
        plot_adjust_flag = 0;% plot replacement process figure 2
        eval = 0; % used to set random seed, using the same seed generate same random environment
        AH
        W
        
        fig_num = 1
        epscur = 1;
        time_vec=[];
        wp_hist = [];

        save_result = false;
        save_video = false;
        videoObj = [];
        video_dt = 0.01;
        video_dt_multiplier = 2;

        t_now = 0;

        plot_fancy_vehicle = false;
        plot_pdf = false;
        pdf_num_sample = 40; % increase this for higher pdf resolutio

        threshold; % risk threshold
        visualize = 1; % turn this off if visualization is not needed
        
        pdf_template = [];
        veh_template = [];
        FRS_plot_data = [];
    end

    methods
        %% constructor
        function S = rlsimulator(AH,W,varargin)
            S = parse_args(S,varargin{:}) ;
            S.AH = AH;
            S.W  = W;
            figure(1);
            set(gcf,'Position',[-491.8000 1.0682e+03 1900 842.4000]); % modify this if needed
            set(gcf,'color','w');
        end
        
        %% step
        function [Observation,Reward,IsDone,LoggedSignals,varargout] = step(S,action)
            %step in time and check for done flags
            tic
            agent_info = S.AH.get_agent_info();
            world_info = S.W.get_world_info(agent_info);%where obstacles are

            
            % move 
            [action_replaced, replace_distance, stuck, k, wp] = S.AH.advanced_move(action,world_info);

            collision = S.check_collision_sim(S.t_now,S.AH.A.time(end),world_info);
            num_static_obstacle = S.W.num_cars - S.W.num_moving_cars-1;
            
            

            if S.visualize
                if ~S.plot_fancy_vehicle
                    S.plot();
                    S.W.plot(S.AH.A.get_agent_info); % plot obstacle vehicles. NOTE: this is not actual obstacle plotting! it's just a place holder    
                    S.AH.plot_A(); % plot ego vehicle. NOTE: this is not actual obstacle plotting! it's just a place holder    
                    scatter(S.AH.wpt_global(1), S.AH.wpt_global(2), 360,'k','x','LineWidth',5);
                else
                    if ~isempty(S.veh_template)
                        S.plot();
                        % generate place holder for vehicles
                        % ego vehicle: 
                        temp = S.veh_template.body;
                        ego.body = plot(temp, 'FaceColor', S.veh_template.ego_body_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                        temp = S.veh_template.roof;
                        ego.roof = plot(temp, 'FaceColor', S.veh_template.ego_body_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                        temp = S.veh_template.windl;
                        ego.windl = plot(temp, 'FaceColor', S.veh_template.ego_wind_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                        temp = S.veh_template.windr;
                        ego.windr = plot(temp, 'FaceColor', S.veh_template.ego_wind_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                        temp = S.veh_template.outline;
                        ego.outline = plot(temp, 'FaceColor', 'none', 'EdgeColor','k', 'FaceAlpha',1,'EdgeAlpha',1);
    
                        % static vehicle:
                        for i = 1:num_static_obstacle
                            x_now = S.W.envCars(i+1,1);
                            y_now = S.W.envCars(i+1,3);
                            temp = S.veh_template.body;
                            temp.Vertices = temp.Vertices + [x_now, y_now];
                            plot(temp, 'FaceColor', S.veh_template.obs_body_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.roof;
                            temp.Vertices = temp.Vertices + [x_now, y_now];
                            plot(temp, 'FaceColor', S.veh_template.obs_body_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.windl;
                            temp.Vertices = temp.Vertices + [x_now, y_now];
                            plot(temp, 'FaceColor', S.veh_template.obs_wind_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.windr;
                            temp.Vertices = temp.Vertices + [x_now, y_now];
                            plot(temp, 'FaceColor', S.veh_template.obs_wind_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.outline;
                            temp.Vertices = temp.Vertices + [x_now, y_now];
                            plot(temp, 'FaceColor', 'none', 'EdgeColor','k', 'FaceAlpha',1,'EdgeAlpha',1);
                        end
    
                        % moving vehicle:
                        for i = 1 : S.W.num_moving_cars
                            temp = S.veh_template.body;
                            temp.Vertices = temp.Vertices - [1000, 1000];
                            obs(i).body = plot(temp, 'FaceColor', S.veh_template.obs_body_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.roof;
                            temp.Vertices = temp.Vertices - [1000, 1000];
                            obs(i).roof = plot(temp, 'FaceColor', S.veh_template.obs_body_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.windl;
                            temp.Vertices = temp.Vertices - [1000, 1000];
                            obs(i).windl = plot(temp, 'FaceColor', S.veh_template.obs_wind_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.windr;
                            temp.Vertices = temp.Vertices - [1000, 1000];
                            obs(i).windr = plot(temp, 'FaceColor', S.veh_template.obs_wind_color, 'EdgeColor','none', 'FaceAlpha',1,'EdgeAlpha',1);
                            temp = S.veh_template.outline;
                            temp.Vertices = temp.Vertices - [1000, 1000];
                            obs(i).outline = plot(temp, 'FaceColor', 'none', 'EdgeColor','k', 'FaceAlpha',1,'EdgeAlpha',1);
                        end
    
                        scatter(S.AH.wpt_global(1), S.AH.wpt_global(2), 360,'k','x','LineWidth',5);
                    end
                end

                cnt = 0;
                while abs(S.t_now - S.AH.A.time(end)) > S.video_dt
                    if S.plot_fancy_vehicle && isempty(S.veh_template)
                        S.plot();
                        scatter(S.AH.wpt_global(1), S.AH.wpt_global(2), 360,'k','x','LineWidth',5);
                    end

                    % get current pos of ego vehicle
                    idx = find(S.AH.A.time <= S.t_now, 1, 'last');
                    alpha = (S.t_now - S.AH.A.time(idx)) / (S.AH.A.time(idx+1) - S.AH.A.time(idx));
                    state_now = S.AH.A.state(:,idx)*(1-alpha) + S.AH.A.state(:,idx+1)*alpha;
                    xy_now = state_now(1:2);
                    h_now = state_now(3);
    
                    % plot pdf
                    if S.plot_pdf
                        mu_sigma_idx = ( 1:1300:(1300*S.W.num_moving_cars-1299) ) + cnt; % this is decided by the data structure of_dynamic_car_world
                        for i = 1:length(mu_sigma_idx)
                            ii = mu_sigma_idx(i);
                            mu = S.W.mu_sigma{1}(ii, 1:2)';
                            Sigma = S.W.mu_sigma{1}(ii, 3:5);
                            Sigma = [(sqrt(Sigma(1))+2.4/3)^2, Sigma(2); Sigma(2), (sqrt(Sigma(3)) + 1.1/3)^2]; % account for footprint ONLY for visualization
                            if isempty(S.pdf_template)
                                if cnt ~= 0
                                    delete(h(i));
                                end
                                pdf = obs_pdf('gaussian','mu',mu,'sigma',Sigma);
                                h(i) = pdf.plot_pdf_environment(S.pdf_num_sample);
                                holder = get(gca,'Children');
                                set(gca,'Children',[holder(2:end-6); holder(1); holder(end-5:end)]);
                            else % NOTICE!!!!! ASSUME Sigma_xy == 0 HERE!!!!! Need to rotate frame if one wants to have nonzero Sigma_xy 
                                xtemp = S.pdf_template.x * sqrt(Sigma(1,1)) + mu(1);
                                ytemp = S.pdf_template.y * sqrt(Sigma(2,2)) + mu(2);
                                if cnt == 0
                                    h(i) = imagesc(xtemp - 100000,ytemp - 10000,S.pdf_template.w,'AlphaData',S.pdf_template.alpha,'AlphaDataMapping','none'); 
                                    colormap(S.pdf_template.colormap);
                                    holder = get(gca,'Children');
                                    set(gca,'Children',[holder(2:end-6); holder(1); holder(end-5:end)]);
                                end
                                if (xy_now(1) - mu(1))^2 <= 250^2 % only plot pdf in the scene to speed up visualiation
                                    h(i).XData = xtemp;
                                    h(i).YData = ytemp;
                                end
                            end
                            
                        end
                        cnt = cnt + S.video_dt_multiplier;
                    end
    
    
                    % ego vehicle
                    if ~S.plot_fancy_vehicle 
                        footprint = [cos(h_now), -sin(h_now); sin(h_now) cos(h_now)]*[-2.4, 2.4, 2.4, -2.4, -2.4; -1.1, -1.1, 1.1, 1.1, -1.1]+xy_now;
                        S.AH.A.plot_data.footprint.Vertices = footprint';
                    else
                        if isempty(S.veh_template)
                            plot_vehicle(xy_now(1), xy_now(2), h_now, [0,0,0]/255, [140,140,140]/255, 1);
                        else
                            rot = [cos(h_now), -sin(h_now); sin(h_now), cos(h_now)];
                            temp = S.veh_template.body;
                            ego.body.Shape.Vertices = temp.Vertices * rot' + [xy_now(1), xy_now(2)];
                            temp = S.veh_template.roof;
                            ego.roof.Shape.Vertices = temp.Vertices * rot' + [xy_now(1), xy_now(2)];
                            temp = S.veh_template.windl;
                            ego.windl.Shape.Vertices = temp.Vertices * rot' + [xy_now(1), xy_now(2)];
                            temp = S.veh_template.windr;
                            ego.windr.Shape.Vertices = temp.Vertices * rot' + [xy_now(1), xy_now(2)];
                            temp = S.veh_template.outline;
                            ego.outline.Shape.Vertices = temp.Vertices * rot' + [xy_now(1), xy_now(2)];
                        end
                    end
                    plot([S.AH.A.state(1,1:idx),xy_now(1)],[S.AH.A.state(2,1:idx),xy_now(2)],'k','LineWidth',2);                   
    
                    % move other vehicles                    
                    if S.plot_fancy_vehicle && isempty(S.veh_template) % plot static vehicle
                        for i = 1:num_static_obstacle
                            x_now = S.W.envCars(i+1,1);
                            y_now = S.W.envCars(i+1,3);
                            if isempty(S.veh_template)
                                plot_vehicle(x_now,y_now,0, [255,255,255]/255, [200,200,200]/255, 1);
                            end
                        end
                    end
                    for i = 1 : S.W.num_moving_cars
                        idx = find(S.W.dyn_obspostime{i} <= S.t_now, 1, 'last');
                        if idx < length(S.W.dyn_obspostime{i})
                            alpha = (S.t_now - S.W.dyn_obspostime{i}(idx)) / (S.W.dyn_obspostime{i}(idx+1) - S.W.dyn_obspostime{i}(idx));
                            x_now = S.W.dyn_obs_pos{i}(idx, 1) * (1-alpha) + S.W.dyn_obs_pos{i}(idx+1, 1) * alpha;
                            y_now = S.W.dyn_obs_pos{i}(idx, 2) * (1-alpha) + S.W.dyn_obs_pos{i}(idx+1, 2) * alpha;
                        else
                            x_now = S.W.dyn_obs_pos{i}(idx, 1);
                            y_now = S.W.dyn_obs_pos{i}(idx, 2);
                        end
                        if S.plot_fancy_vehicle
                            if isempty(S.veh_template)
                                plot_vehicle(x_now,y_now,0, [255,255,255]/255, [200,200,200]/255, 1);
                            else
                                if (x_now - xy_now(1))^2 < 250^2 % only plot obstacles in the scene to speed up visualization
                                    temp = S.veh_template.body;
                                    obs(i).body.Shape.Vertices = temp.Vertices + [x_now, y_now];
                                    temp = S.veh_template.roof;
                                    obs(i).roof.Shape.Vertices = temp.Vertices + [x_now, y_now];
                                    temp = S.veh_template.windl;
                                    obs(i).windl.Shape.Vertices = temp.Vertices + [x_now, y_now];
                                    temp = S.veh_template.windr;
                                    obs(i).windr.Shape.Vertices = temp.Vertices + [x_now, y_now];
                                    temp = S.veh_template.outline;
                                    obs(i).outline.Shape.Vertices = temp.Vertices + [x_now, y_now];
                                end
                            end
                        else
                            xmid = 0.5*( min(S.W.plot_data.obstacles_seen.XData(:,i+num_static_obstacle)) +  max(S.W.plot_data.obstacles_seen.XData(:,i+num_static_obstacle)) );
                            ymid = 0.5*( min(S.W.plot_data.obstacles_seen.YData(:,i+num_static_obstacle)) +  max(S.W.plot_data.obstacles_seen.YData(:,i+num_static_obstacle)) );
                            S.W.plot_data.obstacles_seen.XData(:,i+num_static_obstacle) =  S.W.plot_data.obstacles_seen.XData(:,i+num_static_obstacle) - xmid + x_now;
                            S.W.plot_data.obstacles_seen.YData(:,i+num_static_obstacle) =  S.W.plot_data.obstacles_seen.YData(:,i+num_static_obstacle) - ymid + y_now;
                        end
                    end
    
                    if xy_now(1)+200 <= S.W.goal(1) + 40
                        xlim([xy_now(1)-10, xy_now(1)+200]);
                    else
                        xlim([S.W.goal(1)+40-210, S.W.goal(1)+40]);
                    end

                    
                    ylim([-5, 12]);
                    title("Speed="+num2str(state_now(4),'%.1f')+" [m/s]");
                    S.t_now = S.t_now + S.video_dt_multiplier * S.video_dt;
                    if S.save_video
                        frame = getframe(gcf);
                        writeVideo(S.videoObj, frame);
                    end
                    pause(0.01)
                end
            end


            S.wp_hist{end+1} = wp;
            agent_info = S.AH.get_agent_info() ;
            agent_info.replace_distance = replace_distance;

            agent_info.type = class(S.AH.A);
            Observation =[]; %S.W.get_ob(agent_info);
            Reward = 0;%S.W.getRew(agent_info,Observation);vx = v*cos(h)

            if collision
                warning("A collision Happened!");
            end
            goal_check = S.W.goal_check(agent_info);

            
            IsDone = S.determine_isDone_flag(collision,action_replaced,stuck,goal_check);
            

            
            if S.eval &&( IsDone == 1 || IsDone == 3 ||IsDone == 4 || IsDone == 5) && S.save_result
                Filename = sprintf('Highway_IsDone-%s_SimIdx-%s_SceIdx-%s_thre-%s.mat', num2str(IsDone), num2str(S.W.simulation_idx), num2str(S.epscur), num2str(S.threshold)); 
                envCars = S.W.envCars;
                hist_info.p_hist = S.AH.p_hist;
                hist_info.FRS_hist = S.AH.FRS_hist;
                hist_info.mirror_hist = S.AH.mirror_hist;
                hist_info.type_manu_hist = S.AH.type_manu_hist;
                hist_info.state_hist = S.AH.state_hist;
                hist_info.time_hist = S.AH.time_hist;
                hist_info.wp_hist = S.wp_hist;
                hist_info.solve_time_hist = S.AH.solve_time_hist;                
                traj = S.AH.A.state;
                traj_time = S.AH.A.time;
                save(Filename,'envCars','traj','traj_time','hist_info');
            end
            
            drawnow;
            LoggedSignals = struct;
            
            % send output action if nargout > 4
            if nargout > 4
                varargout = {action_replaced,k} ;
            end
            t_step = toc;
            S.time_vec = [S.time_vec t_step];
        end
        
        function [out] = check_collision_sim(S, tstart, tend, world_info)
            T_check = tstart:0.01:tend;
            out = 0;

            % ego trajectory
            idx1 = max(1,find(S.AH.A.time>=tstart,1));
            ego_time = S.AH.A.time(idx1:end);
            ego_traj = S.AH.A.state(:,idx1:end);

            % static obstacle
            O_all = world_info.obstacles;
            dyn_O.obs   = world_info.dyn_obstacles;
            dyn_O.num_dyn_cars = world_info.num_moving_cars;
            dyn_O.num_static_cars = world_info.num_cars - world_info.num_moving_cars - 1;
            bounds = world_info.bounds; % boundary of the world: [xmin xmax ymin ymax]
            xlo = bounds(1) ; xhi = bounds(2) ;
            ylo = bounds(3) ; yhi = bounds(4) ;
            Blower = [xlo, xhi, xhi, xlo, xlo ; ylo, ylo, ylo-1, ylo-1, ylo] ;
            Bupper = [xlo, xhi, xhi, xlo, xlo ; yhi, yhi, yhi+1, yhi+1, yhi] ;
            B = [Blower, nan(2,1), Bupper, nan(2,1)] ; % make top and bottom bounds into obstacles
            O = [O_all B dyn_O.obs{1}(:,1:6*dyn_O.num_static_cars)] ; % all static obstacles including road boundary 
            num_static = size(O,2)/6;

            % moving obstacle
            obs_traj = {};
            obs_time = {};
            for i = 1:S.W.num_moving_cars
                idx1 = max(1,find(S.W.dyn_obspostime{i}>=tstart,1));
                
                if isempty(idx1) && norm(ego_traj(1:2,1) - S.W.dyn_obs_pos{i}(end,1:2)') >= 250
                    continue
                elseif ~isempty(idx1) && norm(ego_traj(1:2,1) - S.W.dyn_obs_pos{i}(idx1,1:2)') >= 250
                    continue
                end

                if isempty(idx1)  % this means long trajectory for this obstacle is not available because the obstacle is already beyond our goal
                    obs_traj{end+1} = repmat(S.W.dyn_obs_pos{i}(end,1:2)', 1, length(T_check));
                    obs_time{end+1} = T_check;
                else
                    idx2 = find(S.W.dyn_obspostime{i}>tend,1)-1;
                    if isempty(idx2)
                        obs_time{end+1} = S.W.dyn_obspostime{i}(idx1:end);
                        obs_traj{end+1} = S.W.dyn_obs_pos{i}(idx1:end,:)';
                        if S.W.dyn_obspostime{i} < tend
                            obs_time{end} = [obs_time{end}, tend];
                            obs_traj{end} = [obs_traj{end}, obs_traj{end}(:,end)];
                        end
                    else
                        obs_time{end+1} = S.W.dyn_obspostime{i}(idx1:idx2);
                        obs_traj{end+1} = S.W.dyn_obs_pos{i}(idx1:idx2,:)';
                    end
                end
            end

            for t_check = T_check
    
                % ego vehicle location
                h_now = match_trajectories(t_check, ego_time, ego_traj(3,:));
                xy_now = [match_trajectories(t_check, ego_time, ego_traj(1,:)); 
                          match_trajectories(t_check, ego_time, ego_traj(2,:))];
                ego_footprint = [cos(h_now), -sin(h_now); sin(h_now) cos(h_now)]*[-2.4, 2.4, 2.4, -2.4, -2.4; -1.1, -1.1, 1.1, 1.1, -1.1]+xy_now;


                % collision against static obstacle
                for i = 1:num_static
                    obs_footprint = O(:,i*6-5:i*6-1);
                    [xx,yy] = polyxpoly(ego_footprint(1,:),ego_footprint(2,:), obs_footprint(1,:),obs_footprint(2,:));
                    if ~isempty(xx)
                        out = 1;
                        return
                    end

                end

                % collision against moving vehicles
                for i = 1:length(obs_time)
                    obs_traj_check = [match_trajectories(t_check, obs_time{i}, obs_traj{i}(1,:));
                                      match_trajectories(t_check, obs_time{i}, obs_traj{i}(2,:))];
                    obs_footprint = [-2.4, 2.4, 2.4, -2.4, -2.4; -1.1, -1.1, 1.1, 1.1, -1.1]+obs_traj_check;
                    [xx,yy] = polyxpoly(ego_footprint(1,:),ego_footprint(2,:), obs_footprint(1,:),obs_footprint(2,:));
                    if ~isempty(xx)
                        out = 1;
                        return
                    end
                end
            end
        end

        %% reset
        function [iniOb, LoggedSignals] = reset(S)
            LoggedSignals = struct;
            flags = struct;
            flags.discrete_flag = S.discrete_flag;
            flags.replace_action = S.replace_action;
            flags.safety_layer = S.safety_layer;
            if S.eval
                S.W.setup(S.epscur);
                S.AH.reset(flags,S.epscur);
            else
                S.W.setup();
                S.AH.reset(flags);
            end
            iniOb =[];

            if S.save_video
                Videoname = sprintf('Highway_SimIdx-%s_SceIdx-%s_thre-%s', num2str(S.W.simulation_idx), num2str(S.epscur), num2str(S.threshold));
                S.videoObj = VideoWriter(Videoname,'Motion JPEG AVI');
                S.videoObj.FrameRate = 1/(S.video_dt * S.video_dt_multiplier);
                open(S.videoObj);
            end

        end
        
        %% helper functions
        function plot(S)
            if S.plot_sim_flag
                figure(S.fig_num) ;
                cla ; hold on ; axis equal ;

                % plot road
                if ~check_if_plot_is_available(S.W,'road_lanes')
                    road_lanes ={};
                    w= 2;
                    road_lanes{1} =  fill([-S.W.start_line S.W.goal(1)+500 S.W.goal(1)+500 -S.W.start_line -S.W.start_line],[-0.7 -0.7 3*S.W.lanewidth+0.7 3*S.W.lanewidth+0.7 -0.7]-0.5*S.W.lanewidth,[190,190,190]/256); % JL_plot: road 
                    road_lanes{2} =  plot([-S.W.start_line,S.W.goal(1)+500],[3*S.W.lanewidth, 3*S.W.lanewidth]-0.5*S.W.lanewidth,'LineWidth',w,'Color',[255, 255, 255]/255);
                    road_lanes{3} =  plot([-S.W.start_line,S.W.goal(1)+500],[2*S.W.lanewidth, 2*S.W.lanewidth]-0.5*S.W.lanewidth,'--','LineWidth',w,'Color',[1 1 1]);
                    road_lanes{4} =  plot([-S.W.start_line,S.W.goal(1)+500],[S.W.lanewidth, S.W.lanewidth]-0.5*S.W.lanewidth,'--','LineWidth',w,'Color',[1 1 1]);
                    road_lanes{5} =  plot([-S.W.start_line,S.W.goal(1)+500],[0, 0]-0.5*S.W.lanewidth,'LineWidth',w,'Color',[255, 255, 255]/255);
                    S.W.plot_data.road_lanes = road_lanes;
                end

                goal_pos = make_box(S.W.goal_radius*2);
                patch(goal_pos(1,:)+S.W.goal(1),goal_pos(2,:)+S.W.goal(2),[1 0 0]) ;

                % plot FRS
                S.AH.plot_FRS();
                
                xlabel('x [m]');
                ylabel('y [m]');
                set(gca,'FontSize',15)
            end
        end
        
        function [IsDone] = determine_isDone_flag(S,collision,action_replaced,stuck,goal_check)
            if collision && action_replaced
                IsDone = 3;
            elseif collision
                IsDone = 1;
            elseif goal_check
                IsDone = 5;
            elseif stuck
                IsDone = 4;
            elseif action_replaced
                IsDone = 2;
            else
                IsDone = 0;
            end
        end
    end
end
