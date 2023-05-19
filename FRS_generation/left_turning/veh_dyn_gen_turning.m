% states: x, y, h, u, v, r,u0,v0, r0,t0 ,p_u,   p_y,  brake_time, Fy_error, Fx_error r_err_sum, h_err_sum, u_err_sum, h0   t 
%         1  2  3    4  5  6  7  8  9  10  11    12    13         14           15      16        17        18         19    20
% NOTE: states above are adopted from REFINE. HOWEVER, t0(10) and brake_time(13) are no longer in use.

clear, close all, clc
load my_const.mat
syms h brake_time Fy_error Fx_error t p_u u0 t0 v0 r0 u v r err_r_sum err_h_sum err_u_sum h0 
syms x y tdummy udummy p_y Tdummy
dyn = [x; y; h; u; v; r; u0;v0; r0;t0;p_u; p_y;  brake_time; Fy_error; Fx_error; err_r_sum; err_h_sum; err_u_sum; h0; t];


[~,U,Z] = parametrized_turnning_with_brake(p_u, p_y, 1, t); 

%% phase 1: start turning
phase_id = 1;
ud = U(1,phase_id);
dud = U(4,phase_id);
rd = U(3,phase_id);
drd = U(6,phase_id);
hd = Z(phase_id);
ddyn = gen_low_speed(h,u,ud,dud,rd,err_u_sum);
matlabFunction(ddyn, 'File', 'dyn_turn_1lo', 'vars', {tdummy dyn udummy});
ddyn = gen_closed_loop_dyn(h, u, v, r, hd-h0, ud, dud, rd, drd, Fy_error, Fx_error, err_r_sum, err_h_sum, err_u_sum);
matlabFunction(ddyn, 'File', 'dyn_turn_1hi', 'vars', {tdummy dyn udummy});
ddyn = resetmap(dyn,rd);
matlabFunction(ddyn, 'File', 'dyn_turn_1_lo2hi', 'vars', {tdummy dyn udummy Tdummy});

%% phase 2: constant turning
phase_id = 2;
ud = U(1,phase_id);
dud = U(4,phase_id);
rd = U(3,phase_id);
drd = U(6,phase_id);
hd = Z(phase_id);
ddyn = gen_closed_loop_dyn(h, u, v, r, hd-h0, ud, dud, rd, drd, Fy_error, Fx_error, err_r_sum, err_h_sum, err_u_sum);
matlabFunction(ddyn, 'File', 'dyn_turn_2', 'vars', {tdummy dyn udummy});

%% phase 3: back to 0 yaw rate
phase_id = 3;
ud = U(1,phase_id);
dud = U(4,phase_id);
rd = U(3,phase_id);
drd = U(6,phase_id);
hd = Z(phase_id);
ddyn = gen_closed_loop_dyn(h, u, v, r, hd-h0, ud, dud, rd, drd, Fy_error, Fx_error, err_r_sum, err_h_sum, err_u_sum);
matlabFunction(ddyn, 'File', 'dyn_turn_3', 'vars', {tdummy dyn udummy});

%% phase 4: contingency brake before t_stop
phase_id = 4;
ud = U(1,phase_id);
dud = U(4,phase_id);
rd = U(3,phase_id);
drd = U(6,phase_id);
hd = Z(phase_id);
ddyn = gen_closed_loop_dyn(h, u, v, r, hd-h0, ud, dud, rd, drd, Fy_error, Fx_error, err_r_sum, err_h_sum, err_u_sum);
matlabFunction(ddyn, 'File', 'dyn_turn_4', 'vars', {tdummy dyn udummy});

%% phase 5: contingency brake after t_stop
phase_id = 5;
ud = U(1,phase_id);
dud = U(4,phase_id);
rd = U(3,phase_id);
drd = U(6,phase_id);
hd = Z(phase_id);
ddyn = gen_low_speed(h,u,ud,dud,rd,err_u_sum);
matlabFunction(ddyn, 'File', 'dyn_turn_5', 'vars', {tdummy dyn udummy});

%% help functions (dynamics)
function dx = gen_low_speed(h, u, ud, dud, rd, err_u_sum)
    load my_const.mat
    l = lf + lr;
%     Cus = m * grav_const * (lr / (l * Caf1) - lf / (l * Car1));
%     delta = rd / (u+0.05) *(l+Cus*u^2/grav_const); % 0.05 is for numerical issue
    
    mr = lf/l *m;
    
%     rlo = delta*u/(l+Cus*u^2/grav_const);         
    rlo = rd;
    vlo = rlo*(lr - u^2*mr/Car1);
    
    if isSim
        Mu_lo = max_Fx_uncertainty_braking / m;
    else
        Mu_lo = b_u_pro*u+b_u_off;
    end

    kappaU = kappaPU + kappaIU*(err_u_sum);
    phiU = phiPU + phiIU*(err_u_sum);
    u_err = u - ud;
    err_term = Ku*u_err;
    tau_u = -(kappaU*Mu_lo + phiU) * err_term;

    Fxf = - m*vlo*rlo/2 - m/2*Ku*(u-ud) + m/2*dud;
    if isSim
        dudt = 1/m*(2*(Fxf)+m*vlo*rlo) + Mu_lo + tau_u;
    else
        dudt = 1/m*(2*(Fxf)+m*vlo*rlo) + tau_u + (9.8095)*u+(-41.129)*u^2+(60.344)*u^3 ; % (9.8095)*u+(-41.129)*u^2+(60.344)*u^3  is an over approximation of delta_u based on REFINE-Fig.7
    end

    dx = [  u .* cos( h ) - vlo .* sin( h ); u .* sin( h ) + vlo .* cos( h );  rlo; dudt; 0; 0; zeros(9,1);0; 0; 0;0; 1 ];%h_state does not exist

    
    
end
function [dx] = gen_closed_loop_dyn(h, u, v, r, hd, ud, dud, rd, drd, Fy_error, Fx_error, err_r_sum, err_h_sum, err_u_sum)
    load my_const.mat      
    coef = 0.; % something funky happens when Kh is nonzero. Can't tell why. Guess it's a problem of 2pi = 0 ?
    Kh = Kh*coef;

    kappa = kappaP + kappaI*(err_r_sum + 1*err_h_sum) ;
    phi = phiP + phiI*(err_r_sum + 1*err_h_sum);
    r_err = r-rd;
    h_err = h-hd;
    err_term = Kr*r_err + Kh*h_err;          
    tau_r = -(kappa*Mr + phi) *err_term;
    
    kappaU = kappaPU + kappaIU*(err_u_sum);
    phiU = phiPU + phiIU*(err_u_sum);
    u_err = u - ud;
    err_term = Ku*u_err;
    tau_u = -(kappaU*Mu + phiU) * err_term;
    
    dudt = -Ku * ( u - ud )+ dud + Fx_error/m + tau_u;
    drdt = -Kr * ( r - rd )+ drd - Kh * (h - hd) - lr*Fy_error/Izz + tau_r;  % note lf*Fyf-lr*Fyr <= lr*(|Fyf|+|Fyr|) = lr*Fy_error given lr > lf. So here we are giving larger error than it has to be 
    Fyrest =  - Car1 * ( v - lr * r )/(u);
    Fyf =  lr/lf*Fyrest + Izz/lf*( -Kr * ( r - rd )+ drd - Kh * (h - hd) +tau_r);
    dvdt = 1/m * (Fyf +Fyrest + Fy_error) - u*r + Fy_error/m;
    
    dx = [  u .* cos( h ) - v .* sin( h ); u .* sin( h ) + v .* cos( h );  r; dudt; dvdt; drdt; zeros(9,1); r_err^2; h_err^2; u_err^2; 0; 1];
end

function f = resetmap(x,rd)
    u = x(4);
    load my_const.mat
    l = lf+lr;
    rlo = rd;
    mr = lf/l *m;
    vlo = rlo*(lr - u^2*mr/Car1);
    f = x;
    f(5:6) = [vlo;rlo];
end