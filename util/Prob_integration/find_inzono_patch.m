function [in_zono,idx] = find_inzono_patch(zono, xx, yy, xy0, dx, dy, grid_size)
    x0 = xy0(1);
    y0 = xy0(2);
    xy_ctr = [xx(:)'+dx/2; yy(:)'+dy/2];
    
    % transfer patch-buffered zono into its half-space representation
    c = center(zono);
    G = [generators(zono), diag([dx/2, dy/2])];
    B_neg = [-G(2,:); G(1,:)];
    B_pos = (B_neg./vecnorm(B_neg))';
    B = [B_pos; -B_pos];
    b = [B_pos*c + abs(B_pos*G)*ones(size(G,2),1); 
        -B_pos*c + abs(B_pos*G)*ones(size(G,2),1)];

    % check intersection
    idx = find(max(B*xy_ctr-b)<=0); % idx gives the idx of patch
    in_zono = unique(round((xx(idx)-x0)/dx) + round((yy(idx)-y0)/dy)*(grid_size-1));
end