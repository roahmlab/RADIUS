function cost = highway_cost_fun(state,x_des)
            %state = x y h
            %x_des = x y 
            cost = 0.5*sqrt(norm(x_des(1:2) - state(1:2))) + 0.25 * abs(angdiff(x_des(3),state(3)));
%             angle_diff_state = abs(angdiff(x_des(3),state(3))); 
%             cost = norm([x_des(1:2)-state(1:2);angle_diff_state]);
end