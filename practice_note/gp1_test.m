
function [M, S, V] = gp1_test(gpmodel, m, s)
%% Code
% if ~isfield(gpmodel,'induce') || numel(gpmodel.induce)==0
%     [M, S, V] = gp0_kuka_planar_dyn(gpmodel, m, s); return; end

D = length(m); E = size(gpmodel.targets,2);
%% Initialization for robot dynamics
persistent dynamics OPTIONS ctrlfcn u0 par jointlist njoint;
if isempty(dynamics)
    dynamics    = @dynamics_kp_nop;
    OPTIONS     = odeset('RelTol', 1e-2, 'AbsTol', 1e-2);
    ctrlfcn     = str2func('zoh');   
    par.dt = gpmodel.stepsize; par.delay = 0; par.tau = gpmodel.stepsize;
    jointlist   = gpmodel.jointi;
    njoint      = length(jointlist);
    u0          = cell(1,njoint);
end

%% Robot Dynamics
n_span  = gpmodel.n_span;
q       = m(jointlist);
qdot    = m(jointlist + njoint);
tau     = m(end-njoint+1:end);

for j = 1:njoint, u0{j} = @(t)ctrlfcn(tau(j,:),t,par); end
[~, y] = ode45(dynamics, linspace(0,gpmodel.stepsize,n_span+1), m(1:(2*njoint)), OPTIONS, u0{:});

qddot   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,q,qdot,tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
[dqddotdq, dqddotdqdot, dqddotdtau] = solveForwardDynamicsDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,q,qdot,qddot,gpmodel.robot.G,gpmodel.Vdot0,gpmodel.robot.F);
for i = 2:n_span
    q_tmp       = y(i,jointlist)';
    qdot_tmp    = y(i,jointlist + njoint)';
    qddot_tmp   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,q_tmp,qdot_tmp,tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
    [dqddotdq_temp, dqddotdqdot_temp, dqddotdtau_temp] = solveForwardDynamicsDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,q_tmp,qdot_tmp,qddot_tmp,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
    
    dqddotdq            = dqddotdq + dqddotdq_temp;
    dqddotdqdot         = dqddotdqdot + dqddotdqdot_temp;
    dqddotdtau          = dqddotdtau + dqddotdtau_temp;
end
dqddotdq            = dqddotdq * gpmodel.stepsize/n_span;
dqddotdqdot         = dqddotdqdot * gpmodel.stepsize/n_span;
dqddotdtau          = dqddotdtau * gpmodel.stepsize/n_span;

M = zeros(E,1);
M(jointlist)            = (y(n_span+1,jointlist)' - q);
M(jointlist + njoint)   = (y(n_span+1,jointlist + njoint)' - qdot);

A                                        = zeros(E,D);
A(jointlist,jointlist + njoint)          = eye(njoint,njoint) * gpmodel.stepsize;
A(jointlist + njoint,jointlist)          = dqddotdq;
A(jointlist + njoint,jointlist + njoint) = dqddotdqdot;
A(jointlist + njoint,end-njoint+1:end)   = dqddotdtau;




S =  A * s * (A');
S = (S + S')/2;

V = A';
end
function u = zoh(f, t, par) % **************************** zero-order hold
d = par.delay;
if d==0
                  u = f;
else
  e = d/100; t0=t-(d-e/2);
  if t<d-e/2,     u=f(1);
  elseif t<d+e/2, u=(1-t0/e)*f(1) + t0/e*f(2);    % prevents ODE stiffness
  else            u=f(2);
  end
end
end