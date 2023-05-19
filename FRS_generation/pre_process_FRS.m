function res = pre_process_FRS(FRS, manu_type, grid_size)
    % this function transfers a zonotope into the data structure that Cuda 
    % can use for probability integration
    % FRS: a cell vector of zonotopes
    % manu_type: maneuver type. 1(spd), 2(dir), 3(lan)
    % grid_size: grid to cover FRS is expected to be grid_size X grid_size.
    %            default value 25

    if nargin < 3
        grid_size = 25;
    end

    %% CUDA information
    res.grid_size = grid_size - 1;
    res.num_zono = length(FRS);
    res.grid_x0 = [];
    res.grid_y0 = [];
    res.grid_dx = [];
    res.grid_dy = [];
    res.g_u0_x = []; % x dimension of the generators corresponding to u0
    res.g_u0_y = []; % y dimension of the generators corresponding to u0
    res.g_v0_x = []; % x dimension of the generators corresponding to v0
    res.g_v0_y = []; % y dimension of the generators corresponding to v0
    res.g_r0_x = []; % x dimension of the generators corresponding to r0
    res.g_r0_y = []; % y dimension of the generators corresponding to r0
    res.block_inzono_list = [];
    res.rot_angle = [];
    res.g_p_x = []; % xy dimensions of the generator corresponding to p
    res.cg_p = []; % center and generator on the dimension of p


    %% c++ information for Lulu
    res.manu_type = manu_type;
    res.cg_u0v0r0 = []; % center and generator on the dimensions of u0 v0 r0 
    res.t_range = [];

    %% collect data
    for i = 1:length(FRS)

        bla = pre_process_zono(FRS{i},manu_type,grid_size);
         
        res.grid_x0 = [res.grid_x0, bla.grid_xy0(1)]; % [x0_zono1, x0_zono2, ...]
        res.grid_y0 = [res.grid_y0, bla.grid_xy0(2)]; % [y0_zono1, y0_zono2, ...]
        res.grid_dx = [res.grid_dx, bla.grid_dx]; % [dx_zono1, dx_zono2, ...]
        res.grid_dy = [res.grid_dy, bla.grid_dy]; % [dy_zono1, dy_zono2, ...]
        res.g_u0_x = [res.g_u0_x, bla.g_u0v0r0(1,1)]; % [g_u0_x_zono1, g_u0_x_zono2, ...]
        res.g_u0_y = [res.g_u0_y, bla.g_u0v0r0(2,1)]; % [g_u0_y_zono1, g_u0_y_zono2, ...] 
        res.g_v0_x = [res.g_v0_x, bla.g_u0v0r0(1,2)]; % [g_v0_x_zono1, g_v0_x_zono2, ...]
        res.g_v0_y = [res.g_v0_y, bla.g_u0v0r0(2,2)]; % [g_v0_y_zono1, g_v0_y_zono2, ...]
        res.g_r0_x = [res.g_r0_x, bla.g_u0v0r0(1,3)]; % [g_r0_x_zono1, g_r0_x_zono2, ...]
        res.g_r0_y = [res.g_r0_y, bla.g_u0v0r0(2,3)]; % [g_r0_y_zono1, g_r0_y_zono2, ...]
        res.rot_angle = [res.rot_angle, bla.rot_angle]; % [rot_angle_zono1, rot_angle_zono2, ...]
        res.g_p_x = [res.g_p_x, bla.g_p(1)]; % [g_p_x_zono1, g_p_x_zono2, ...]
        res.t_range = [res.t_range, bla.t_range']; % [tmin_zono1, tmax_zono1, tmin_zono2, tmax_zono2]

        % block_inzono_list is a bit different. We use grid_size^2 blocks
        % to cover each zonotope. For the block that intersects with the
        % block, we save its index. Otherwise we assign -1. Thus we save 
        % (grid_size-1)^2 numbers for zonotope coverage. However, because
        % in Cuda each zonotope is computed by blocks of blockDim.x=32
        % and blockDim.y = blockDim.z = 1.Therefore for coaliscing, we may 
        % save n numbers for each zonotope where mod(n,32) = 0 and 
        % n >= grid_size^2.
        n = ceil((grid_size - 1)^2 / 32) * 32;
        num_neg1 = n - length(bla.block_inzono_list);
        res.block_inzono_list = [res.block_inzono_list, bla.block_inzono_list, -ones(1,num_neg1)]; 
        % [block_inzono_list_zono1, block_inzono_list_zono2, ...].
        % NOTE: length(block_inzono_list_zono1) = length(block_inzono_list_zono2) = n, a CONSTANT!!!
    end

    res.cg_p = [res.cg_p, bla.c_p, bla.g_p(3)]; % [c_p_zono1, g_p_zono1]. All zonotopes in FRS share the same cg_p 
    res.cg_u0v0r0 = [res.cg_u0v0r0, bla.c_u0v0r0', bla.g_u0v0r0(end,:)]; % [c_u0v0r0_zono1 (1-by-3), g_u0v0r0_zono1 (1-by-3)]. All zonotopes in FRS share the same info 
 
end