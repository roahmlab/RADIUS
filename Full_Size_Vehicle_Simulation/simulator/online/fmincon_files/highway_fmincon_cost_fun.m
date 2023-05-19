function [cost,dcost] = highway_fmincon_cost_fun(K,agent_state,slice_idx,FRS,x_des)
% %K,agent_state, slice_idx,FRS.vehRS_save{FRS.brake_idx(1)}{1},x_des
% [zono_slice, betas, slice_generators,slice_G_full]=zonotope_slice_cost_fun(FRS, [7;8;9;slice_idx], [agent_state(4);agent_state(5);agent_state(6);K]);
% future_pose = center(zono_slice);
% cost = highway_cost_fun(future_pose(1:3),x_des);
% delta = future_pose(1:2)-x_des(1:2); 
% cost = sum(delta.^2);
% dcost = sum(2.*(delta).*slice_G_full(1:2,end)/slice_generators(4,4));


% JL cost function
idx = find(FRS.Z(slice_idx,2:end))+1;
slice_cg = FRS.Z([1,2,slice_idx],[1,idx]);
lambda = (K-slice_cg(3,1))/slice_cg(3,2);
delta = slice_cg(1:2,1) + slice_cg(1:2,2)*lambda - x_des(1:2);
cost = [1,50] * delta.^2; % JL
dcost = 2*[1,50] * ( delta.*slice_cg(1:2,2)/slice_cg(3,2) );
% beta_1 = 
% dcost = 0;

if isnan(cost) || isnan(dcost)
    cost
end

end
