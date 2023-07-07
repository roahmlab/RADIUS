function [T U Z] = parametrized_turnning_with_brake(p_u, p_y, symbolic_flag, t)
% we assume we will start with a very small initial speed condition, thus we will keep accelerating during the first phase with amax no matter what

    load my_const.mat

    % The speed parameter value should not be less than the critical speed,
    % this most likely indicates an input error
    assert(u_cri <= p_u);

    tm1 = tpk_dir/3;
    tm2 = 2;
    tm3 = tm1; % tm = tm1+tm2+tm3
    tb1 = max([(p_u-u_cri)/amax, 0]);
    tb2 = 1;
    dt = 0.001;
    if symbolic_flag
        t1 = t; % start turning
        t2 = t; % constant turning
        t3 = t; % yaw rate goes to 0
        t4 = t; % contigency braking before t_stop
        t5 = t; % contigency braking after t_stop
        T = [];
    else
        t1 = 0:dt:tm1;
        t2 = (dt:dt:tm2)+tm1;
        t3 = (dt:dt:tm3)+tm1+tm2;
        t4 = (dt:dt:tb1)+tm1+tm2+tm3;
        t5 = (dt:dt:tb2)+tm1+tm2+tm3+tb1;
        T = [t1 t2 t3 t4 t5];
    end

    %% start turning
    ud1 = (amax+0.5) * t1;
    dud1 = (amax+0.5) * ones(size(t1));
    vd1 = zeros(size(t1)); % place holder
    dvd1 = vd1;
    rd1  = p_y/2 - (p_y*cos((2*pi*t1)/(tm1+tm3)))/2;
    drd1 = (p_y*pi*sin((2*pi*t1)/(tm1+tm3)))/(tm1+tm3);
    %% constant turning
    ud2 = ones(size(t2)) * p_u;
    dud2 = zeros(size(t2));
    vd2 = zeros(size(t2)); 
    dvd2 = vd2;
    rd2 = p_y/2 - (p_y*cos((2*pi*tm1*ones(size(t2)))/(tm1+tm3)))/2;
    drd2 = vd2;
    %% yaw rate goes to 0
    ud3 = ones(size(t3)) * p_u;
    dud3 = zeros(size(t3));
    vd3 = zeros(size(t3)); 
    dvd3 = vd3;
    rd3 = p_y/2 - (p_y*cos((2*pi*(t3 - tm2))/(tm1+tm3)))/2;
    drd3 = (p_y*pi*sin((2*pi*(t3-tm2))/(tm1+tm3)))/(tm1+tm3);
    %% contigency braking before t_stop
    ud4 = p_u - amax * (t4-tm1-tm2-tm3);
    dud4 = -amax * ones(size(t4));
    vd4 = zeros(size(t4));
    dvd4 = vd4;
    rd4 = vd4;
    drd4 = vd4;
    %% contigency braking after t_stop
    ud5 = zeros(size(t5));
    dud5 = ud5;
    vd5 = ud5;
    dvd5 = vd5;
    rd5 = vd5;
    drd5 = vd5;

    %% Group
    ud = [ud1 ud2 ud3 ud4 ud5];
    dud = [dud1 dud2 dud3 dud4 dud5];
    vd = [vd1 vd2 vd3 vd4 vd5];
    dvd = [dvd1 dvd2 dvd3 dvd4 dvd5];
    rd = [rd1 rd2 rd3 rd4 rd5];
    drd = [drd1 drd2 drd3 drd4 drd5];
    U = [ud; vd; rd; dud; dvd; drd];
    if ~symbolic_flag
        h = cumsum(rd)*dt;
        x_dot = ud.*cos(h) - vd.*sin(h);
        y_dot = ud.*sin(h) + vd.*cos(h);
        Z = [cumsum([x_dot;y_dot],2)*dt; h];
    else
        h1 = int(rd1,t1);
        h2 = subs(int(rd2,t2),t2,t2-tm1) + subs(h1,t1,tm1);
        h3 = int(rd3,t3) + rd2*tm2 - rd2;
        h4 = p_y*(tm1+tm3)/2 + rd2*tm2;
        h5 = h4;
        Z = [h1 h2 h3 h4 h5];
%         [-(p_y*(sin(pi*t) - pi*t))/(2*pi), (p_y*(2*t - 1))/2, (p_y*(2*pi - sin(pi*t) + pi*t))/(2*pi), 3*p_y, 3*p_y]
    end
    

end


