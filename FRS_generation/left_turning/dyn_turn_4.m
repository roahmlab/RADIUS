function ddyn = dyn_turn_4(tdummy,in2,udummy)
%DYN_TURN_4
%    DDYN = DYN_TURN_4(TDUMMY,IN2,UDUMMY)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-May-2023 15:18:51

Fx_error = in2(15,:);
Fy_error = in2(14,:);
err_h_sum = in2(17,:);
err_r_sum = in2(16,:);
err_u_sum = in2(18,:);
h = in2(3,:);
h0 = in2(19,:);
p_u = in2(11,:);
p_y = in2(12,:);
r = in2(6,:);
t = in2(20,:);
u = in2(4,:);
v = in2(5,:);
t2 = cos(h);
t3 = sin(h);
t4 = p_u.*4.0;
t5 = u.*4.0;
t6 = t.*2.0e+1;
t9 = err_h_sum.*1.005102352581729;
t10 = err_r_sum.*1.005102352581729;
t7 = -t5;
t8 = -t6;
t11 = t9+t10+4.002551176290865;
mt1 = [t2.*u-t3.*v;t3.*u+t2.*v;r;Fx_error./1.575e+3+t4+t7+t8+(err_u_sum.*(3.28e+2./3.75e+2)+1.624380952380952).*(t4+t7+t8+8.0e+1)+7.5e+1;Fy_error.*1.26984126984127e-3-r.*3.678044669195112+((r.*4.847676e+5-v.*2.9028e+5).*1.573254670599803e-3)./u-r.*t11.*3.678044669195112-r.*u];
mt2 = [Fy_error.*(-5.1023525817293e-4)-r.*2.0-r.*t11.*2.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;r.^2;(h+h0-p_y.*3.0).^2;(p_u-t.*5.0-u+2.0e+1).^2;0.0;1.0];
ddyn = [mt1;mt2];
end
