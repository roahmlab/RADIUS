function res = pre_process_zono(ZONO,manu_type,grid_size)
    % this function transfers a zonotope into the data structure that Cuda 
    % can use for probability integration
    % ZONO: a sliceable zonotope generated from CORA with dimension:
    %           [x,y,h,u,v,r,u0,v0,r0,...,Pu,Py,...,t]
    %            1,2,3,4,5,6, 7, 8, 9,...,11,12,...,end
    % manu_type: maneuver type. 1(spd), 2(dir), 3(lan)
    % grid_size: grid to cover FRS is expected to be grid_size X grid_size.
    
    if nargin < 3
        grid_size = 25;
    end

    Z = ZONO.Z;
    dim = size(Z,1);
    if manu_type == 1
        P_dim = 11;
    else
        P_dim = 12;
    end

    %%
    % buffer FRS by other vehicle's footprint. Assuming other vehicle
    % travels straight with speed upto 22 m/s. Vehicle size = 2.2*4.8
    % ATTENTION: CHANGE L IF WE KNOW OUR AND MOVING VEHICLE'S HEADINGS!!!!!
    t_itv = interval(zonotope(Z(dim,:)));
    t_range = [infimum(t_itv); supremum(t_itv)];
%     L = norm([2.2, 4.8]) + 25*(supremum(t_itv) - infimum(t_itv)); % for hardware, need to change vehicle size!!!
%     Z = [Z, [L/2*eye(2); zeros(dim-2,2)]];
    
%     hitv = interval(zonotope(Z(3,:))); % for simulation
%     hrad = rad(hitv) + 0.05*0;
%     hmid = -mid(hitv);


    % aggressive bounding
%     if supremum(t_itv) - infimum(t_itv) == 0.01
%         L = 4.8 + 25*0.01;
%     else
%         L = 4.8 + 25*0.02;
%     end
    L = 4.8 + 0 * (supremum(t_itv) - infimum(t_itv)); % account displacement during the time interval in pdf % turns out we don't neet to!
    W = 2.2;
    hrad = 0.03*0;

%     SL = sqrt(L^2 + W^2);
%     SW = L*abs(sin(hrad)) + W*abs(cos(hrad));
%     ROT = [cos(hmid), -sin(hmid); sin(hmid) cos(hmid)];
    
    SL = [cos(hrad), -sin(hrad)] * [L;-W];
    SW = [sin(hrad),  cos(hrad)] * [L; W];
    ROT = eye(2);

    Z = [Z, [ROT*diag([-SL/2, -SW/2]); zeros(dim-2,2)]];

%     plot(zonotope(Z),[1,2],'g');

    %%
    % deal with sliceable dimension (7,8,9 and P_dim), and kick these
    % sliceable generator out
    c_slice = Z([7,8,9,P_dim],1);
    g_slice = [];
    for i = [7,8,9,P_dim]
        idx = find(Z(i,2:end)~=0)+1;
        g_slice = [g_slice, Z([1,2,i],idx)]; % [x;y;sliceable_dimension]
        Z(:,idx) = [];
        if isempty(idx)
            % in case no sliceable generator is found, add fake sliceable
            % but tiny generator. this is because old generator doesn't
            % have initial error on v0 and r0
            g_slice = [g_slice,ones(3,1)*1e-10];
        end
    end

    % rotate the zonotope w.r.t. the sliceable generator w.r.t. P. 
    % PS: pdf needs to rotate the same angle as well
    th = atan2(g_slice(2,4), g_slice(1,4)); 
    R = [cos(th), -sin(th); sin(th), cos(th)]';
    Z = R*Z(1:2,:);
    g_slice(1:2,:) = R*g_slice(1:2,:);

    % cover the zonotope by grid_size X grid_size grid, and check
    % intersection of each small block
    zono = zonotope(Z);
    zono_itv = interval(zono);
    zono_bound = [infimum(zono_itv), supremum(zono_itv)];
    [xx,yy] = meshgrid(linspace(zono_bound(1,1),zono_bound(1,2),grid_size), linspace(zono_bound(2,1),zono_bound(2,2),grid_size));
    grid_dx = xx(1,2) - xx(1,1);
    grid_dy = yy(2,1) - yy(1,1);
    xy0 = zono_bound(:,1);
    
    [in_zono, idx_inzono] = find_inzono_patch(zono,xx(1:grid_size-1,1:grid_size-1),yy(1:grid_size-1,1:grid_size-1),xy0,grid_dx,grid_dy,grid_size);
%%%%%%% below is a naive brute-force implementation of the line above %%%%%%%%%%    
%     in_zono = [];
%     for i = 1:grid_size-1
%         for j = 1:grid_size - 1
%             small_itv = interval([xx(i,j);yy(i,j)], [xx(i,j);yy(i,j)]+[grid_dx;grid_dy]);
%             if isIntersecting(zono,small_itv)
%                 in_zono = [in_zono, (i-1)*(grid_size-1)+(j-1)]; % index of the block
%                 % blocks in the grid are index as:
%                 % |       ...,         ..., ...,           ...|
%                 % | grid_size, grid_size+1, ..., 2*grid_size-1|
%                 % |         0,           1, ...,   grid_size-1|
%                 continue
%             end
% 
%             %%% CORA's isIntersecting function is wrong because it has
%             %%% false negative !!! Add sanity check here
%             pts = [xx(i,j), xx(i,j) + grid_dx, xx(i,j) + grid_dx, xx(i,j);
%                    yy(i,j), yy(i,j), yy(i,j) + grid_dy, yy(i,j) + grid_dy];
%             for k = 1:4
%                 if containsPoint(zono, pts(:,k))
%                     in_zono = [in_zono, (i-1)*(grid_size-1)+(j-1)];
%                     break
%                 end
%             end
%         end
%     end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % save result
    res.c_xy = Z(1:2,1);           % 2 by 1
    res.c_u0v0r0 = c_slice(1:3);   % 3 by 1
    res.g_u0v0r0 = g_slice(:,1:3); % 3 by 3
    res.c_p = c_slice(4);          % 1 by 1
    res.g_p = g_slice(:,4);        % 3 by 1
    res.t_range = t_range;         % 2 by 1
    res.grid_xy0 = xy0;            % 2 by 1
    res.grid_dx = grid_dx;         % 1 by 1
    res.grid_dy = grid_dy;         % 1 by 1
    res.block_inzono_list = in_zono; % 1 by n
    res.rot_angle = -th;           % 1 by 1
    res.rotated_buf_zono = zono;


end