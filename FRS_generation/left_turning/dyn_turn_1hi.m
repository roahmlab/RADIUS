function ddyn = dyn_turn_1hi(tdummy,in2,udummy)
%DYN_TURN_1HI
%    DDYN = DYN_TURN_1HI(TDUMMY,IN2,UDUMMY)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-May-2023 15:18:48

Fx_error = in2(15,:);
Fy_error = in2(14,:);
err_h_sum = in2(17,:);
err_r_sum = in2(16,:);
err_u_sum = in2(18,:);
h = in2(3,:);
h0 = in2(19,:);
p_y = in2(12,:);
r = in2(6,:);
t = in2(20,:);
u = in2(4,:);
v = in2(5,:);
t2 = cos(h);
t3 = sin(h);
t4 = t.*pi;
t5 = r.*2.0;
t6 = u.*4.0;
t9 = -p_y;
t10 = t.*2.2e+1;
t13 = err_h_sum.*1.005102352581729;
t14 = err_r_sum.*1.005102352581729;
t7 = cos(t4);
t8 = sin(t4);
t11 = -t6;
t16 = t13+t14+4.002551176290865;
t12 = p_y.*t7;
t15 = t5+t9+t12;
et1 = Fy_error.*1.26984126984127e-3+p_y.*1.839022334597556-r.*3.678044669195112-t12.*1.839022334597556+((r.*4.847676e+5-v.*2.9028e+5).*1.573254670599803e-3)./u-r.*u-t15.*t16.*1.839022334597556;
et2 = p_y.*t8.*pi.*9.19511167298778e-1;
ddyn = [t2.*u-t3.*v;t3.*u+t2.*v;r;Fx_error./1.575e+3+t10+t11-(err_u_sum.*(3.28e+2./3.75e+2)+1.624380952380952).*(t6-t10)+1.1e+1./2.0;et1+et2;Fy_error.*(-5.1023525817293e-4)+p_y-t5+t7.*t9-t15.*t16+(p_y.*t8.*pi)./2.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;(p_y.*(-1.0./2.0)+r+t12./2.0).^2;(h+h0-(p_y.*(t4-t8))./(pi.*2.0)).^2;(t.*(1.1e+1./2.0)-u).^2;0.0;1.0];
end