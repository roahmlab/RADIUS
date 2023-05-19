function d = min_dist_between_ego_and_dyn_obs(...
    ego_min_v, ...
    ego_max_v, ...
    obs_v, ...
    ego_x0, ...
    obs_x0, ...
    t_range_sec)
    delta_x_init = (ego_x0 - obs_x0);
    delta_vel_min = (ego_min_v - obs_v);
    delta_vel_max = (ego_max_v - obs_v);
    get_min_dist_time = @(delta_vel, delta_x) -(delta_x ./ delta_vel);
    get_dist_at_time = @(t, delta_vel, delta_x) abs(t .* (delta_vel) + delta_x);
    clamp_t = @(t) arrayfun(@(x) max([min([x, t_range_sec]), 0]), t)
    t_min_vel_min = clamp_t(get_min_dist_time(delta_vel_min, delta_x_init));
    t_min_vel_max = clamp_t(get_min_dist_time(delta_vel_max, delta_x_init));
    d_min_vel_min = get_dist_at_time(t_min_vel_min, delta_vel_min, delta_x_init);
    d_min_vel_max = get_dist_at_time(t_min_vel_max, delta_vel_max, delta_x_init);
    d = min([d_min_vel_min], [d_min_vel_max]);
end