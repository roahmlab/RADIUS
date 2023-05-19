function [T U Z] = sin_one_hump_parameterized_traj_with_brake(t0,del_y,p_u,u0,t,symbolic_flag, not_speed_change,brake_time)
% this function provides desired trajectory for a speed or direction change maneuver
load my_const.mat
if(not_speed_change)
     p_y = del_y;
else
     p_y = 0; 
end

if symbolic_flag
    t1=t;
    t2=t;
    t3=t;
    T = [];
    Z = [];
    tm = tpk_dir;
    if ~exist('brake_time')
        t_stop = (p_u - u_really_slow)/amax; 
    else
        t_stop = brake_time;
    end
    
else 
    t_stop = (p_u - u_really_slow)/amax; 
    if t_stop < 0
        t_stop = 0;
    end
    tm = tpk_dir;
    tf = tm + t_stop + tbrk2 -t0; 
    T = linspace(0, tf,1000);
    Z = [];
    idx = find(T>=(tm-t0),1)-1;
    t1 = T(1:idx);
    idx2 = find(T>=(tm+t_stop-t0),1)-1;
    t2 = T((idx+1):idx2);
    t2 = t2 - (tm-t0);
    t3 = T(idx2+1:end);
    t3 = t3 - (tm+t_stop-t0);
end

%% specify ud, vd , rd during the driving maneuver portion
t1_shifted = t1 + t0;
ud1 = (p_u-u0)/tm*t1_shifted+u0; 
dud1 = (p_u-u0)/tm*ones(1,length(t1_shifted));

vd1 = zeros(size(t1_shifted));
dvd1 =zeros(size(t1_shifted));
rd1  = p_y/2 - (p_y*cos((2*pi*t1_shifted)/tm))/2;
drd1 = (p_y*pi*sin((2*pi*t1_shifted)/tm))/tm;

%% contingency braking before t_stop
% (Do not use the symbolic mode to generate braking trajectory!!!!!)
if symbolic_flag % place holder
    ud2 = p_u - (p_u-u_really_slow)*(t2)/t_stop;
    dud2 = -(p_u-u_really_slow)/t_stop*ones(1,length(t2));
    vd2 = zeros(1,length(t2));
    dvd2 = zeros(1,length(t2));
    rd2 = zeros(1,length(t2));
    drd2 = zeros(1,length(t2));
else
    if t_stop == 0
        ud2 = []; dud2 = []; vd2 = []; dvd2 = []; rd2 = []; drd2 =[];
    else
        ud2 = p_u - (p_u-u_really_slow)*(t2)/t_stop;
        dud2 = -(p_u-u_really_slow)/t_stop*ones(1,length(t2));
        vd2 = zeros(1,length(t2));
        dvd2 = zeros(1,length(t2));
        rd2 = zeros(1,length(t2));
        drd2 = zeros(1,length(t2));
    end
end
%% contingency braking after t_stop
t3_shifted = t3;
ud3 = 0*t3_shifted;
dud3 =0*t3_shifted;

vd3 = zeros(1,length(t3));
dvd3 = zeros(1,length(t3));
rd3 = zeros(1,length(t3));
drd3 = zeros(1,length(t3));




%% Group
ud = [ud1 ud2 ud3];
dud=[dud1 dud2 dud3];
vd = [vd1 vd2 vd3];
dvd = [dvd1 dvd2 dvd3];
rd  = [rd1 rd2 rd3];
drd  = [drd1 drd2 drd3];
U = [ud;vd;rd;dud;dvd;drd];
if ~symbolic_flag
    %integrate desired as reference and output as Z
    p = cumsum([rd1 rd2 rd3])*(t1(2)-t1(1));
    x_dot = ud.*cos(p) - vd.*sin(p);
    y_dot = ud.*sin(p) + vd.*cos(p);
    Z = cumsum([x_dot;y_dot],2)*(t1(2)-t1(1));
    Z = [Z;p];
else
    Z = [int(rd1,t)-subs(int(rd1,t),t, 0);0;0];
end
end