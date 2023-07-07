function mu_sigmas = dyn_obs_to_mu_sigmas(x0, y0, h0, ...
	vel, len_local_wrt_obs, width_local_wrt_obs, ...
	dt_seconds, t_width_seconds, t_total_seconds, ...
	lanewidth, t_now, ...
	sigma_x_local_fcn, sigma_y_local_fcn)
    % Generates a set of (mu_x, mu_y, sigma_xx, sigma_xy, sigma_yy) values over a given set of time intervals
    % 
    % INPUTS:
    %  x0: the x coordinate center of the obstacle at t = 0
    %  y0: the y coordinate center of the obstacle at t = 0
    %  h0: the h coordinate center of the obstacle at t = 0
    %  vel: the velocity of the obstacle
    %  len_local_wrt_obs: the length of the obstacle (thickness along its local x-axis)
    %  width_local_wrt_obs: the width of the obstacle (thickness along its local y-axis)
    %  dt_seconds: the time delta between each mu sigma
    %  t_width_seconds: the duration that each mu sigma covers
    %  lanewidth: the width of a lane
    %  t_now: the start time, offset from t = 0
    %  sigma_x_local_fcn: (local_length, local_width, lanewidth) -> positive real number
    %  sigma_y_local_fcn: (local_length, local_width, lanewidth) -> positive real number
    num_out_dbl = t_total_seconds / dt_seconds - ((t_width_seconds / dt_seconds) - 1);
    num_out = int64(num_out_dbl);
    assert(abs(num_out_dbl - double(num_out)) < 1.0e-19, ...
        "Must have integer number of mu sigmas")
    mu_sigmas = [];
    start_center_t = t_width_seconds / 2;
    new_length = t_width_seconds * vel + len_local_wrt_obs;
    rotmat_h0 = rotmat(h0);
    for i = 1:num_out
        t_center = double(i-1) * dt_seconds + start_center_t + t_now;
        d_center_from_start = vel * t_center;
        mu_x = x0 + cos(h0) * d_center_from_start;
        mu_y = y0 + sin(h0) * d_center_from_start;
        %sigma_x_local = (new_length/2/3)^2;
        %sigma_y_local = ((1.32)/2/3)^2;
        sigma_x_local = sigma_x_local_fcn(new_length, width_local_wrt_obs, lanewidth);
        sigma_xy_local = 0;
        sigma_y_local = sigma_y_local_fcn(new_length, width_local_wrt_obs, lanewidth);
        sigma_mat_local = [...
            sigma_x_local  sigma_xy_local; ...
            sigma_xy_local sigma_y_local   ...
        ];
        sigma_mat_global = rotmat_h0 * sigma_mat_local * rotmat_h0';
        
        mu_sigmas(end+1:end+5) = [ ...
            mu_x,                  ...
            mu_y,                  ...
            sigma_mat_global(1,1), ...
            sigma_mat_global(1,2), ...
            sigma_mat_global(2,2)  ...
        ];
    end
end
