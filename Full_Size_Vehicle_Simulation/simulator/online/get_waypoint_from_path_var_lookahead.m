function [wp,wp1,wp2] = get_waypoint_from_path_var_lookahead(vehicle_pos, path, max_lookahead)
    [d,~,~,p_idx] = dist_point_on_polyline_simon(vehicle_pos(1:2),path(1:2,:)) ;
    d_best = dist_polyline_cumulative(path(1:2,:));
%     reduce_curve_lookahead = 7;
    if (p_idx + 10)< length(path(1,:))%20 points = 1m 
        vert = path(1:2,p_idx:p_idx+10)';
        k=LineCurvature2D(vert)/20;
        N=LineNormals2D(vert);
        figure(1);
        plot([vert(:,1) vert(:,1)+k.*N(:,1)]',[vert(:,2) vert(:,2)+k.*N(:,2)]','m','LineWidth',3);
%         plot([vert(Lines(:,1),1) vert(Lines(:,2),1)]',[vert(Lines(:,1),2) Vertices(Lines(:,2),2)]','b');
        k = max(abs(LineCurvature2D(path(1:2,p_idx:p_idx+10)')));
        if k > 0.2
            lookahead_distance = 3*15;% 1.5 s peak velocity, 0.55 m/s
%         lookahead_distance = lookahead_distance_potential+min_lookahead;
        else
            lookahead_distance = 3*15;% 1.5 s peak velocity, 2m/s
        end 
    else
        lookahead_distance = 3*6;
% %         "max lookahead"
% %         lookahead_distance = 1;
    end
    
%     lookahead_distance =1.5*0.55
    
   
    d_lkhd = min(d_best(end), d + lookahead_distance);
    wp = match_trajectories(d_lkhd,d_best,path);
    d_lkhd1 = min(d_best(end), d + lookahead_distance/3);
    wp1 = match_trajectories(d_lkhd1,d_best,path);
    d_lkhd2 = min(d_best(end), d + lookahead_distance/3*2);
    wp2 = match_trajectories(d_lkhd2,d_best,path);
    
    
    %if you want heading lookahead to be a little further to get through
    %corners, then use
%     heading_lookahead = lookahead_distance + 0.8;
%     d_lkhd = min(d_best(end), d + heading_lookahead);
%     wp_h = match_trajectories(d_lkhd,d_best,path);
%     d_lkhd1 = min(d_best(end), d + heading_lookahead/3);
%     wp_h1 = match_trajectories(d_lkhd1,d_best,path);
%     d_lkhd2 = min(d_best(end), d + heading_lookahead/3*2);
%     wp_h2 = match_trajectories(d_lkhd2,d_best,path);
%     
%     wp(3) = wp_h(3);
%     wp1(3) = wp_h1(3);
%     wp2(3) = wp_h2(3);
end