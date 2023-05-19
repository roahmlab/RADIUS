function zono = combine_zono(zono1, zono2, manu_type, time_flag)
    % manu_type = 1 for spd, 2 for dir, 3 for lan

    if nargin < 4
        time_flag = 0;
    end

    if manu_type == 1
        special_idx = [7 8 9 11];
    elseif manu_type == 4
        special_idx = 12;
    else
        special_idx = [7 8 9 12];
    end

    Z1 = zono1.Z; 
    Z2 = zono2.Z;

    % deal with special generators first
    G1_special = [];
    G2_special = [];
    for i = special_idx
        idx = find(Z1(i,2:end))+1;
        G1_special = [G1_special, Z1(:,idx)];
        Z1(:,idx) = [];

        idx = find(Z2(i,2:end))+1;
        G2_special = [G2_special, Z2(:,idx)];
        Z2(:,idx) = [];
    end
    G_special = 0.5*(G1_special + G2_special);
    G_special_add = 0.5*(G1_special - G2_special);

    c1 = Z1(:,1); c2 = Z2(:,1);
    c = 0.5*(c1+c2); 
    c_add = 0.5*(c1-c2);
    
    % get time information (not generalizd!!!)
    tmid = c1(20);
    if time_flag
        if abs(round(tmid*100) - tmid*100) > 1e-4
            t1_range = [tmid-0.005, tmid+0.005];
        else
            t1_range = [tmid-0.01, tmid+0.01];
        end
    else
        bla = interval(zono1);
        bla = bla(20);
        t1_range = [infimum(bla), supremum(bla)];
    end
    tmid = c2(20);
    if time_flag
        if abs(round(tmid*100) - tmid*100) > 1e-4
            t2_range = [tmid-0.005, tmid+0.005];
        else
            t2_range = [tmid-0.01, tmid+0.01];
        end
    else
        bla = interval(zono2);
        bla = bla(20);
        t2_range = [infimum(bla), supremum(bla)];
    end
    t_range = [min([t1_range, t2_range]), max([t1_range, t2_range])];
    


    G1 = Z1(:,2:end); G2 = Z2(:,2:end);

    % uniform sign
    G1 = uniform_sign(G1);
    G2 = uniform_sign(G2);

    % prune small generators
    G1 = prune_gen(G1); G2 = prune_gen(G2);
    idx = find(G1(1,:).^2 + G1(2,:).^2 <= 0.0001^2);
    G1_too_small = G1(:,idx); G1(:,idx) = [];
    idx = find(G2(1,:).^2 + G2(2,:).^2 <= 0.0001^2);
    G2_too_small = G2(:,idx); G2(:,idx) = [];


    % match generators
    G1 = [G1; atan2( round(G1(2,:)*1000)/1000, round(G1(1,:)*1000)/1000 ); vecnorm(G1(1:2,:))]; % add information for angles
    G2 = [G2; atan2( round(G2(2,:)*1000)/1000, round(G2(1,:)*1000)/1000 ); vecnorm(G2(1:2,:))];
    G1 = sortrows(G1', [22,21], 'descend')'; G2 = sortrows(G2', [22,21], 'descend')';
%     idx = find(G1(21,:) == 0);
%     if ~isempty(idx)
%         G1(:,idx) = sortrows(G1(:,idx)',1,'descend')';
%     end
%     idx = find(G2(21,:) == 0);
%     if ~isempty(idx)
%         G2(:,idx) = sortrows(G2(:,idx)',1,'descend')';
%     end
%     idx = find(G1(21,:) == pi/2);
%     if ~isempty(idx)
%         G1(:,idx) = sortrows(G1(:,idx)',2,'descend')';
%     end
%     idx = find(G2(21,:) == pi/2);
%     if ~isempty(idx)
%         G2(:,idx) = sortrows(G2(:,idx)',2,'descend')';
%     end

    G1_match = []; G1_unmatch = []; 
    G2_match = []; G2_unmatch = G2; 
    for i = 1:size(G1,2)
        a1 = G1(21,i);
%         [val, idx] = min(abs(G2_unmatch(21,:) - a1));
        if min(abs(G2_unmatch(21,:) - a1)) > 0.1
            G1_unmatch = [G1_unmatch, G1(:,i)];
        else
            idx = find(abs(G2_unmatch(21,:) - a1) <= 0.1);
            [~,bla] = min(abs(G2_unmatch(22,idx) - G1(22,i)));
            idx = idx(bla);
            G1_match = [G1_match, G1(:,i)];
            G2_match = [G2_match, G2_unmatch(:,idx)];
            G2_unmatch(:,idx) = [];
        end
    end


    % combine all generators
    if isempty(G1_match) && isempty(G2_match)
        G = []; G_add = [];
    else
        if size(G1_match,2) <= size(G2_match,2)
            G = 0.5*(G1_match(1:20,:) + G2_match(1:20,1:size(G1_match,2)));
            G_add = [0.5*(G1_match(1:20,:) - G2_match(1:20,1:size(G1_match,2))), G2_match(1:20,size(G1_match,2)+1:end)];
        else
            G = 0.5*(G2_match(1:20,:) + G1_match(1:20,1:size(G2_match,2)));
            G_add = [0.5*(G2_match(1:20,:) - G1_match(1:20,1:size(G2_match,2))), G1_match(1:20,size(G2_match,2)+1:end)];
        end
    end
    Z = [c, G_special, G, c_add, G_special_add, G_add];


    if ~isempty(G1_unmatch)
        Z = [Z, G1_unmatch(1:20,:)];
    end
    if ~isempty(G2_unmatch)
        Z = [Z, G2_unmatch(1:20,:)];
    end

    Z(20,:) = Z(20,:) * 0;
    Z(20,1) = mean(t_range);
    Z(20,1+size(G_special,2)+1) = 0.5*(t_range(2) - t_range(1));


    zono = deleteZeros(deleteAligned_noslice(zonotope(Z), [special_idx,20]));
end





function M = uniform_sign(M)
    % unify the sign of first two rows
    idx = find(sign(M(1,:)) == -1);
    M(1:2,idx) = -M(1:2,idx); % want the sign of first dim to be non-negative

    idx = find(sign(M(1,:)) == 0);
    idx2 = find(sign(M(2,:)) == -1);
    idx = intersect(idx, idx2);
    M(1:2,idx) = -M(1:2,idx); % when the first dim is 0, want the sign of second dim to be non-negative


    % deal with 0 dimensions
    idx = find(abs(M(1,:)) <= 1e-4);
    if ~isempty(idx)
        M(:,idx) = M(:,idx) .* my_sign(M(2,idx));
        bla = M(:,1)*0;
        bla(1) = sum(abs(M(1,idx)));
        bla(2) = sum(M(2,idx));
        M(:,idx) = [];
        M = [M, bla];
    end

    idx = find(abs(M(2,:)) <= 1e-4);
    if ~isempty(idx)
        M(:,idx) = M(:,idx) .* my_sign(M(1,idx));
        bla = M(:,1)*0;
        bla(2) = sum(abs(M(2,idx)));
        bla(1) = sum(M(1,idx));
        M(:,idx) = [];
        M = [M, bla];
    end
end

function a = my_sign(a)
    if a == 0
        a = 1;
    else
        a = sign(a);
    end
end

function M = prune_gen(M)
    idx = find(M(1,:).^2 + M(2,:).^2 <= 0.001^2);
    m = M(:,idx);
    M(:,idx) = [];
    M = [M, sum(abs(m),2)];
end