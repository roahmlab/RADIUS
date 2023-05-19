clear, clc;

u0 = 0.1;
p_u = 10;
p_y = 0.52;
delpy = 0.1;

load my_const.mat
dim = 20;
options.tensorParallel = 0;%turn this to 1 if want to use parallel not working very well.
options.taylorTerms=15; % number of taylor terms for reachable sets
options.zonotopeOrder= 100; % zonotope order... increase this for more complicated systems.
options.maxError = 1e100*ones(dim, 1); 
options.verbose = 1;
options.uTrans = 0; % we won't be using any inputs
options.U = zonotope([0, 0]);
options.advancedLinErrorComp = 0;
options.tensorOrder = 2;
options.errorOrder = 50;
options.reductionInterval = inf;
options.reductionTechnique = 'girard';
options.tStart = 0;
options.timeStep = 0.01;
options.tStart = 0;
options.tFinal = 0;
options.x0 = zeros(dim,1);
options.x0([4,7]) = u0; % initial speed
options.x0(11) = p_u; %p_u
options.x0(12) = p_y;
delp_y = zeros(dim,1); delp_y(12) = delpy;
delFy = zeros(dim,1); delFy(14) = max_Fy_uncertainty_braking;
delFx =  zeros(dim,1);  delFx(15) = max_Fx_uncertainty_braking;

vehRs = [];
%% phase 1
tnow = 0;
tlength = 0.01;
Rt{1}{1}.set = zonotope([options.x0, delp_y, delFy, delFx]);
t_guard1 = 0;
t_guard2 = 0;
while tnow < 1
    bla = Rt{1}{1}.set;
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_1lo, Rt{1}{1}.set, options, tlength);

    % check if guard is triggered because in phase 1 we speed up from 0.1 to 5.5
    u = interval(Ri{1}{1});
    u = u(4);
    if ~t_guard1 && supremum(u) > u_cri
        t_guard1 = tnow;
        R_before_guard = bla;
    end
    if infimum(u) >= u_cri
        t_guard2 = tnow; % so we know guard is triggered within [t_guard1, t_guard2]
        break
    end
    vehRs = [vehRs, Ri{1}];
    tnow = tnow + tlength
end

% pass the guard
assert(t_guard1 > 0);
if abs(t_guard2 - t_guard1 - 0.01) < 1e-5 % 1e-5 due to numerical issue. 
    R_lo = vehRs{floor(t_guard2 / tlength)}; % this is all possible states that trigger the guard
else
    options.tFinal = floor(t_guard1 / tlength) * tlength;
    Rt{1}{1}.set = R_before_guard;
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_1lo, Rt{1}{1}.set, options, floor((t_guard2-t_guard1)/tlength)*tlength, floor((t_guard2-t_guard1)/tlength)*tlength);
    R_lo = Ri{1}{1}; % this is all possible states that trigger the guard
end
tnow = t_guard1;
options.tFinal = tnow;
options.timeStep = tlength;

% reset map 
options_aux.zonotopeOrder=100;
options_aux.reachabilitySteps=1; 
options_aux.tensorOrder = 2;
options_aux.errorOrder = 50;
options_aux.reductionTechnique='girard';
options_aux.uTrans = 0;
options_aux.U = zonotope([0,0]);
options_aux.tStart = 0;
options_aux.tFinal = 1;
options_aux.timeStep = 1;
options_aux.uTrans = 0;
options_aux.R0 = R_lo;
options_aux.x0 = center(R_lo);
sys_aux = nonlinearSysDT(dim,1,@dyn_turn_1_lo2hi,options_aux);
R_hi = reach(sys_aux,options_aux);
Rt{1}{1}.set = deleteAligned(deleteZeros(R_hi{2}));


% keep going till the end of phase 1.
while tnow < 1 - options.timeStep
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_1hi, Rt{1}{1}.set, options, tlength);
    vehRs = [vehRs, {Rt{1}{1}.set}];
    tnow = tnow + tlength
end

%% phase 2 and 3
% NOTICE: because we don't know the exact time when the guard is triggered,
% we sve Rt in the previous while loop instead of Ri. However vehRs{end}
% should contain all possible states at time t=1, thus we concentrate time
% and returns to use Ri.
Z = Rt{1}{1}.set.Z;
Z(end,:) = [1, zeros(1,size(Z,2)-1)];

% increase force uncertainty
idx = find(Z(14,:)); Z(14,idx) = max_Fy_uncertainty;
idx = find(Z(15,:)); Z(15,idx) = max_Fx_uncertainty/3;

Rt{1}{1}.set = zonotope(Z);

while tnow < 3 - options.timeStep
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_2, Rt{1}{1}.set, options, tlength);
    vehRs = [vehRs, {Ri{1}{1}.set}];
    tnow = tnow + tlength
end

while tnow < 4 - options.timeStep
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_3, Rt{1}{1}.set, options, tlength);
    vehRs = [vehRs, {Ri{1}{1}.set}];
    tnow = tnow + tlength
end

%% phase 4
tbrake = (p_u - u_cri)/amax + 0.5;
t_guard1 = 0;
t_guard2 = 0;

% decrease force uncertainty
Z = Rt{1}{1}.set.Z;
idx = find(Z(14,:)); Z(14,idx) = max_Fy_uncertainty_braking;
idx = find(Z(15,:)); Z(15,idx) = max_Fx_uncertainty_braking;
Rt{1}{1}.set = zonotope(Z);

while tnow < 4 + tbrake - options.timeStep
    bla = Rt{1}{1}.set;
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_4, Rt{1}{1}.set, options, tlength);
    % check if guard is triggered because in phase 4 we slow down from p_u
    u = interval(Ri{1}{1});
    u = u(4);
    if ~t_guard1 && infimum(u) < u_cri
        t_guard1 = tnow;
        R_before_guard = bla;
    end
    if supremum(u) <= u_cri
        t_guard2 = tnow; % so we know guard is triggered within [t_guard1, t_guard2]
        break
    end
    vehRs = [vehRs, Ri{1}];
    tnow = tnow + tlength
end

% pass the guard
assert(t_guard1 > 0);
if abs(t_guard2 - t_guard1 - 0.01) < 1e-5 % 1e-5 due to numerical issue. 
    R_lo = vehRs{end}; % this is all possible states that trigger the guard
else
    options.tFinal = floor(t_guard1 / tlength) * tlength;
    Rt{1}{1}.set = R_before_guard;
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_4, Rt{1}{1}.set, options, floor((t_guard2-t_guard1)/tlength)*tlength, floor((t_guard2-t_guard1)/tlength)*tlength);
    R_lo = Ri{1}{1}; % this is all possible states that trigger the guard
end
tnow = t_guard1;
options.tFinal = tnow;
options.timeStep = tlength;

%% phase 5 
% back to use Rt again due to the guard issue
Rt{1}{1}.set = R_lo;
while tnow < 5.5 - options.timeStep
    [Ri, Rt, options] = reach_local(dim, @dyn_turn_5, Rt{1}{1}.set, options, tlength);
    vehRs = [vehRs, {Rt{1}{1}.set}];
    tnow = tnow + tlength
end



%% add necessary info for CUDA computation
LeftTurn_FRS2CUDA;





function [Ri, Rt, options] = reach_local(dim, dyn, R0, options, tlength, timeStep)
% states: x, y, h,  u, v, r,u0,v0, r0,t0 ,p_u,  p_y,  brake_time, Fy_error, Fx_error r_err_sum, h_err_sum, u_err_sum, h0   t 
%         1  2  3    4  5  6  7  8  9  10  11    12    13         14           15      16        17        18         19    20
% NOTE: states above are adopted from REFINE. HOWEVER, t0(10) and brake_time(13) are no longer in use.
    c = R0.Z(:,1);
    h0 = c(3);
    xy0 = c(1:2);
    options.R0 = ROTxy(R0, -h0);
    options.x0 = center(options.R0);
    options.tStart = options.tFinal;
    options.tFinal = options.tStart + options.timeStep;
    if exist('timeStep','var')
        options.timeStep = timeStep;
    end
    if exist('tlength','var')
        options.tFinal = options.tStart + tlength;
    end
    options.x0(19) = h0;
    sys = nonlinearSys(dim,1,dyn,options);
    [Ri,Rt] = reach(sys,options); % Ri and Rt have the same length

    % rotate everything back
    for i = 1:length(Ri)
        Ri{i}{1} = deleteAligned(deleteZeros(ROTxy(Ri{i}{1}, h0, xy0)));
        Rt{i}{1}.set = deleteAligned(deleteZeros(ROTxy(Rt{i}{1}.set, h0, xy0)));
    end

end


