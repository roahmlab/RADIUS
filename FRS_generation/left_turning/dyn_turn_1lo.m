function ddyn = dyn_turn_1lo(tdummy,in2,udummy)
%DYN_TURN_1LO
%    DDYN = DYN_TURN_1LO(TDUMMY,IN2,UDUMMY)

%    This function was generated by the Symbolic Math Toolbox version 9.3.
%    19-May-2023 15:18:47

err_u_sum = in2(18,:);
h = in2(3,:);
p_y = in2(12,:);
t = in2(20,:);
u = in2(4,:);
t2 = cos(h);
t3 = sin(h);
t4 = t.*pi;
t5 = u.*4.0;
t6 = u.^2;
t8 = t.*2.2e+1;
t10 = p_y./2.0;
t7 = cos(t4);
t9 = -t5;
t13 = t6.*2.189696155436131e-3;
t12 = p_y.*t7.*(-1.0./2.0);
t15 = t13-1.67e+2./1.0e+2;
t14 = t10+t12;
ddyn = [t2.*u+t3.*t14.*t15;t3.*u-t2.*t14.*t15;t14;t8+t9-(err_u_sum.*(7.0./1.0e+1)+1.3e+1./1.0e+1).*(t5-t8)+1.1e+1./2.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;0.0;1.0];
end
