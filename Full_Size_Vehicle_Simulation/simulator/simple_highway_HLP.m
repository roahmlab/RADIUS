classdef simple_highway_HLP < high_level_planner
% HLP stands for 'high-level planner'
    properties
        lookahead = 60; % lookahead distance
        lanewidth = 4;
    end
methods
    function HLP = simple_highway_HLP(varargin)
        HLP@high_level_planner(varargin{:}) ;
    end
    function reset(HLP)
        HLP.waypoints = [];
        HLP.current_waypoint_index = 1;
    end
    
    function [waypoint,desired_lane] = get_waypoint(HLP,world_info,agent_state)
        % major function for getting the waypoint during planning

%         O_dyn = world_info.dyn_obstacles{1};
%         %reconstruct cars;
%         car_info = zeros(2,0);
%         for i = 1:size(O_dyn,2)/6
%             idx = (i-1)*6+1:i*6-2;
%             car_info = [car_info [min(O_dyn(1,idx));mean(O_dyn(2,idx));world_info.dyn_obstacles{2}(i)]]; % x,y,spd
%         end

        % JL add to replace lines above and let it work for both refine and risk-rtd
        car_info = [world_info.envCars(2:end,1)' + world_info.envCars(2:end,2)' * world_info.time;
                    world_info.envCars(2:end,3)';
                    world_info.envCars(2:end,2)'];
   
        
        car_info = [car_info [99999;0;0] [99999;HLP.lanewidth;0] [99999;HLP.lanewidth*2;0]];
        car_relative_dist = car_info(1,:) - agent_state(1); 
        car_info = [car_relative_dist(:,car_relative_dist>0); car_info(2:3,car_relative_dist>0)];

        car_info(1,:) = car_info(1,:)+6*car_info(3,:) - 10; % JL: consider future position
        car_info(1,car_info(3,:) == 0) = car_info(1,car_info(3,:) == 0) - 120; % JL penalize on static obstacles
        
        car_info_lane_1_logi = car_info(2,:) <= 1;
        car_info_lane_2_logi = car_info(2,:) <= 6 & car_info(2,:) >= 1;
        car_info_lane_3_logi = car_info(2,:) >= 6;
        [lane_1_score]=min(car_info(1,car_info_lane_1_logi));
        [lane_2_score]=min(car_info(1,car_info_lane_2_logi));
        [lane_3_score]=min(car_info(1,car_info_lane_3_logi));
        lane_scores = [lane_1_score lane_2_score lane_3_score];
        [~,desired_lane]=max(lane_scores);
        ego_lane = abs(agent_state(2)-[0,HLP.lanewidth,2*HLP.lanewidth]) < 2.5;
        if any(lane_scores(ego_lane) < 160) && (agent_state(1)>100) % decress the lookahead distance if the ahead vehicle is close
            lookahead_d = 42; 
        else
            lookahead_d = HLP.lookahead; 
        end

        waypoint = [agent_state(1)+lookahead_d; HLP.lanewidth*(desired_lane-1); 0];
%         plot(waypoint(1), waypoint(2),'k*')

        % update current waypoints
        HLP.current_waypoint = waypoint ;
        HLP.waypoints = [HLP.waypoints, waypoint] ;
        HLP.current_waypoint_index = HLP.current_waypoint_index + 1 ;
    end
end
end