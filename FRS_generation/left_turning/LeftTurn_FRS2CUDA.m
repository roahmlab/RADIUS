round_num = 2;
manu_type = 4; % 4 = left turning
grid_size = 33;

frs = vehRs;
figure, hold on 
for i = 1:length(frs)
    plot(frs{i},[1,2],'b');
end

% enclose zonotopes to reduce the number of zonotopes
for j = 1:round_num
    frs_temp = {};
    cnt = 1;
    while cnt < length(frs)
        z = combine_zono(frs{cnt}, frs{cnt+1}, 4);
        intz = interval(z); intz = intz(1:2);
        range_z = supremum(intz) - infimum(intz);
        int1 = interval(frs{cnt}); int1 = int1(1:2);
        int2 = interval(frs{cnt+1}); int2 = int2(1:2);
        mini = min([infimum(int1),infimum(int2)],[],2);
        maxi = max([supremum(int1),supremum(int2)],[],2);
        range_ori = maxi - mini;
        if all(range_z./range_ori <= 1.05)
            frs_temp{end+1} = z;
            cnt = cnt + 2;
        else
            frs_temp{end+1} = frs{cnt};
            cnt = cnt + 1;
        end
    end
    frs = frs_temp;
end
L = 4.8;
W = 2.2;
for i = 1:length(frs)
    plot(frs{i},[1,2],'g');
    hitv = interval(frs{i});
    hitv = hitv(3);
    hmid = mid(hitv);
    hrad = rad(hitv);
    
    rot = [cos(hrad), -sin(hrad); sin(hrad) cos(hrad)];
    pt = rot * [L; -W]/2;

    if pt(2) > 0
        S = diag([sqrt(L^2+W^2), L*abs(sin(hrad)) + W*abs(cos(hrad))])/2;
    else
        S = diag([pt(1), L*abs(sin(hrad))/2 + W*abs(cos(hrad))/2]);
    end
    S = [cos(hmid), -sin(hmid); sin(hmid), cos(hmid)] * S;
    frs{i} = zonotope([frs{i}.Z, [S;zeros(18,2)], [zeros(6,3);0.001*eye(3);zeros(11,3)]]); % add footprint and redundant (auxillary) sliceable generator
    plot(frs{i},[1,2],'m');
end

LeftTurnFRS(1).vehRS_save = frs;
LeftTurnFRS(1).manu_type = 4;
LeftTurnFRS(1).cuda_FRS = pre_process_FRS(frs, manu_type, grid_size);
LeftTurnFRS(1).u_final = 10;

save CUDA_LeftTurn_frs LeftTurnFRS

