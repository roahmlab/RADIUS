function world_obs = get_world_bounds_as_dyn_obs(world_x_min, world_x_max, world_y_min, world_y_max)
% 
% world_x_min = -10;
% world_x_max = 1000;
% 
% world_y_min = 0;
% world_y_max = 10;

world_x_range = world_x_max - world_x_min;
world_y_range = world_y_max - world_y_min;

world_y_mean = world_y_range/2 + world_y_min;
world_x_mean = world_x_range/2 + world_x_min;

world_bound_eps = 0.05;
left_bound_obs = [world_x_min - world_bound_eps; ... % X0
                  world_y_mean; ...                  % Y0
                  0; ...                             % Heading
                  0; ...                             % Vel
                  world_bound_eps * 2; ...           % Length
                  world_y_range * 2; ...             % Width
                 ];
right_bound_obs = [world_x_max + world_bound_eps; ... % X0
                   world_y_mean; ...                  % Y0
                   0; ...                             % Heading
                   0; ...                             % Vel
                   world_bound_eps * 2; ...           % Length
                   world_y_range * 2; ...             % Width
                  ];
top_bound_obs = [world_x_mean; ...                    % X0
                 world_y_max + world_bound_eps; ...   % Y0
                 0; ...                               % Heading
                 0; ...                               % Vel
                 world_x_range * 2; ...               % Length
                 world_bound_eps * 2; ...             % Width
                ];
bottom_bound_obs = [world_x_mean; ...                    % X0
                    world_y_min - world_bound_eps; ...   % Y0
                    0; ...                               % Heading
                    0; ...                               % Vel
                    world_x_range * 2; ...               % Length
                    world_bound_eps * 2; ...             % Width
                   ];
world_obs = [left_bound_obs right_bound_obs ...
             top_bound_obs bottom_bound_obs];
end