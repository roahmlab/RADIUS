classdef highway_HLP < high_level_planner
    properties
        mode = 2; % mode 1 avoid min, mode 2 chase max, other for fixed point
        desired_y = 6;
        T = 3;
    end
methods
    function HLP = highway_HLP(varargin)
        HLP@high_level_planner(varargin{:}) ;
    end
    function reset(HLP)
        HLP.waypoints = [];
        HLP.current_waypoint_index = 1;
    end
    
    function waypoint = get_waypoint(HLP,world_info,agent_state)
        % d is a 3x1 matrix that holds the front car info.
        O_dyn = world_info.dyn_obstacles{1};
        %reconstruct cars;
        car_info = zeros(2,0);
        for i = 1:size(O_dyn,2)/6
            idx = [(i-1)*6+1:i*6-2];
            car_info = [car_info [min(O_dyn(1,idx));mean(O_dyn(2,idx));world_info.dyn_obstacles{2}(i)]];
        end
        car_info = [car_info [99999;0;100] [99999;3.7;100] [99999;3.7*2;100]];
        car_info(1,:) = car_info(1,:) - agent_state(1);
        car_info = car_info(:, car_info(1,:)>0);
        time_head =  car_info(1,:) ./(agent_state(4) - car_info(3,:));
        time_head(time_head<0) = 1000;
        car_info_lane_1_logi = car_info(2,:) <= 1;
        car_info_lane_2_logi = car_info(2,:) <= 6 & car_info(2,:) >= 1;
        car_info_lane_3_logi = car_info(2,:) >= 6;
        
        lan1_cars = time_head(car_info_lane_1_logi);
        lan2_cars = time_head(car_info_lane_2_logi);
        lan3_cars = time_head(car_info_lane_3_logi);
        [lane_1_score,lan1_car_idx]=min(lan1_cars(1,:));
        [lane_2_score,lan2_car_idx]=min(lan2_cars(1,:));
        [lane_3_score,lan3_car_idx]=min(lan3_cars(1,:));
        
%         if lan1_cars(3,lan1_car_idx) < 0.1
%             lane_1_score = 0;
%         end
%         if lan2_cars(3,lan2_car_idx) < 0.1
%             lane_2_score = 0;
%         end
%         if lan3_cars(3,lan3_car_idx) < 0.1
%             lane_3_score = 0;
%         end
        lane_scores = [lane_1_score lane_2_score lane_3_score];
        [~,desired_lane]= max(lane_scores);
        [~, ego_lane]   = min(abs(agent_state(2)-[0,3.7,7.4]));
        if lane_scores(ego_lane) < 10
            lookahead_d = 35; %50; %30;
        else
            lookahead_d = 55; %90; %50;
        end
        waypoint = [agent_state(1)+lookahead_d; 3.7*(desired_lane-1); 0];

        % update current waypoints
        HLP.current_waypoint = waypoint ;
        HLP.waypoints = [HLP.waypoints, waypoint] ;
        HLP.current_waypoint_index = HLP.current_waypoint_index + 1 ;

    end
end
end