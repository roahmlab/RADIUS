function [ddir] = rough_manu_cost_estimate(tbdir,agent_state, x_des,t0_idx)
tbdir = tbdir(:,:,t0_idx);
ddir = [];
homo_robot = homo_mat(agent_state(3),agent_state(1:2));

for i = 1: size(tbdir,2) %Since Ay is one sided, do opposite side with mirror
    homo_manuver = homo_mat(tbdir(4,i),tbdir(2:3,i));
    homo_after = homo_robot*homo_manuver;
    if size(tbdir,1) > 8
    homo_manuver1 = homo_mat(tbdir(11,i),tbdir(9:10,i));
    homo_after1 = homo_robot*homo_manuver1;
    homo_manuver2 = homo_mat(tbdir(14,i),tbdir(12:13,i));
    homo_after2 = homo_robot*homo_manuver2;
    
%     scatter(homo_after(1,3),homo_after(2,3),'b','Filled');
%     axis equal
    rough_cost = highway_cost_fun([homo_after(1:2,3);wrapToPi(tbdir(4,i)+agent_state(3))],x_des(:,1)) + ...
                 highway_cost_fun([homo_after1(1:2,3);wrapToPi(tbdir(11,i)+agent_state(3))],x_des(:,2))+...
                  highway_cost_fun([homo_after2(1:2,3);wrapToPi(tbdir(14,i)+agent_state(3))],x_des(:,3));
    else
        rough_cost = highway_cost_fun([homo_after(1:2,3);wrapToPi(tbdir(4,i)+agent_state(3))],x_des(:,1));
    end
              
    ddir       = [ddir rough_cost]; % heading difference from 0 in rad
end