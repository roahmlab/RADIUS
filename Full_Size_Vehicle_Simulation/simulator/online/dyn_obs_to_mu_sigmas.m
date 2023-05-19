function mu_sigmas = dyn_obs_to_mu_sigmas(x0, y0, h0, vel, len, width, dt_seconds, t_width_seconds, t_total_seconds, lanewidth)
    % TODO should adjust length based on velocity
    assert(h0 == 0, "Currently everything else will break if h0 != 0");
    num_out_dbl = t_total_seconds / dt_seconds - ((t_width_seconds / dt_seconds) - 1);
    num_out = int64(num_out_dbl);
    assert(abs(num_out_dbl - double(num_out)) < 1.0e-19, ...
        "Must have integer number of mu sigmas")
    mu_sigmas = [];
    start_center_t = t_width_seconds / 2;
    new_length = t_width_seconds * vel + len;
    for i = 1:num_out
        t_center = double(i-1) * dt_seconds + start_center_t;
        d_center_from_start = vel * t_center;
        mu_x = x0 + cos(h0) * d_center_from_start;
        mu_y = y0 + sin(h0) * d_center_from_start;
        % sqrt(sigma_x)*3 = length
        % sigma_x = (length/3)^2
        sigma_x = (new_length/2/3)^2;
        sigma_xy = 0;
%         sigma_y = (width/3)^2;
        % Bounding box is width is given by 3 * std
        % = 3 * sqrt( ((lanewidth - 2.2) / 2 / 3)^2 )
        % = 3 * ((lanewidth - 2.2) / 2 / 3)
        % = (lanewidth - 2.2) / 2
        % = (lanewidth - 2.2) / 2
        % TODO 2.2 = width
%         assert(width == 2.2)
%         assert(lanewidth == 3.7)
        sigma_y = ((1.32)/2/3)^2;
        mu_sigmas(end+1:end+5) = [mu_x, mu_y, sigma_x, sigma_xy, sigma_y];
    end
end