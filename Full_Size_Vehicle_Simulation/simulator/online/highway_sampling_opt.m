function [K]= highway_sampling_opt(FRS,O, agent_state, mirror_flag, type_manu)            
FRS = FRS{diridx};
[safeK] = AH.find_safe_actionset(O,agent_state,FRS,mirror_flag,type_manu);
if ~isempty(safeK)
    cost_arr = zeros(length(safeK),1);
    gensize = (tbdir(1,1)-tbdir(1,2))/2;
    safeK = safeK.*gensize+tbdir(1,diridx);
    for i = 1:length(safeK)
        c = center(zonotope_slice(FRS.vehRS_save{FRS.brake_idx(1)}{1}, [7;8;9;12], [agent_state(4);agent_state(5);agent_state(6);safeK(i)]));
        cost_arr(i) = highway_cost_fun(c(1:3),x_des);
    end
    [~,Kidx] = min(cost_arr);
    K = [agent_state(4);safeK(Kidx)];
    break;
else
    K = [];
    if mirror_flag == 0
        ddir(diridx) = 1000;
    else
        ddir_mirror(diridx) = 1000;
    end
end