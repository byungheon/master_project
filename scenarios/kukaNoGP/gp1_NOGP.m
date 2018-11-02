%% gp1_kuka_planar_dyn.m
% *Summary:* Compute joint predictions for the FITC sparse approximation to 
% multiple GPs with uncertain inputs. 
% Predictive variances contain uncertainty about the function, but no noise. 
% If gpmodel.nigp exists, individual noise contributions are added.
%
%   function [M, S, V] = gp1d(gpmodel, m, s)
% 
% *Input arguments:*
%
%   gpmodel    GP model struct
%     hyp      log-hyper-parameters                                  [D+2 x  E ]
%     inputs   training inputs                                       [ n  x  D ]
%     targets  training targets                                      [ n  x  E ]
%     nigp     (optional) individual noise variance terms            [ n  x  E ]
%   m          mean of the test distribution                         [ D  x  1 ]
%   s          covariance matrix of the test distribution            [ D  x  D ]
%
% *Output arguments:*
%
%   M          mean of pred. distribution                            [ E  x  1 ]
%   S          covariance of the pred. distribution                  [ E  x  E ]
%   V          inv(s) times covariance between input and output      [ D  x  E ]
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-05
%
%% High-Level Steps
% # If necessary, compute kernel matrix and cache it
% # Compute predicted mean and inv(s) times input-output covariance
% # Compute predictive covariance matrix, non-central moments
% # Centralize moments

function [M, S, V] = gp1_NOGP(gpmodel, m, s)
%% Code
% if ~isfield(gpmodel,'induce') || numel(gpmodel.induce)==0
%     [M, S, V] = gp0_kuka_planar_dyn(gpmodel, m, s); return; end

D = length(m); E = size(gpmodel.targets,2);
%% Initialization for robot dynamics
persistent dynamics OPTIONS ctrlfcn u0 par jointlist njoint ncontrol;
if isempty(dynamics)
    dynamics    = @dynamics_NOGP_learning;
    OPTIONS     = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
    ctrlfcn     = str2func('zoh');   
    par.dt = gpmodel.stepsize; par.delay = 0; par.tau = gpmodel.stepsize;
    jointlist   = gpmodel.jointi;
    njoint      = length(jointlist);
    ncontrol    = njoint - 1;
    u0          = cell(1,ncontrol);
end


%% Robot Dynamics
n_span  = gpmodel.n_span;
q       = m(jointlist);
qdot    = m(jointlist + njoint);
tau     = [m(end-ncontrol+1:end);0];

for j = 1:ncontrol, u0{j} = @(t)ctrlfcn(tau(j,:),t,par); end
[~, y] = ode45(dynamics, linspace(0,gpmodel.stepsize,n_span+1), m(1:(2*njoint)), OPTIONS, u0{:});

% 
qddot   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q),stateTojoint(qdot),tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
[dqddotdq, dqddotdqdot, dqddotdtau] = solveForwardDynamicsDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q),stateTojoint(qdot),qddot,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
for i = 2:n_span
    q_tmp       = y(i,jointlist)';
    qdot_tmp    = y(i,jointlist + njoint)';
    qddot_tmp   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q_tmp),stateTojoint(qdot_tmp),tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
    [dqddotdq_temp, dqddotdqdot_temp, dqddotdtau_temp] = solveForwardDynamicsDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q_tmp),stateTojoint(qdot_tmp),qddot_tmp,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
    
    dqddotdq            = dqddotdq + dqddotdq_temp;
    dqddotdqdot         = dqddotdqdot + dqddotdqdot_temp;
    dqddotdtau          = dqddotdtau + dqddotdtau_temp;
end
dqddotdq            = dqddotdq * gpmodel.stepsize/n_span;
dqddotdqdot         = dqddotdqdot * gpmodel.stepsize/n_span;
dqddotdtau          = dqddotdtau(:,1:ncontrol) * gpmodel.stepsize/n_span;

M = zeros(E,1);
M(jointlist)            = (y(n_span+1,jointlist)' - q);
M(jointlist + njoint)   = (y(n_span+1,jointlist + njoint)' - qdot);

A                                        = zeros(E,D);
A(jointlist,jointlist + njoint)          = eye(njoint,njoint) * gpmodel.stepsize;
A(jointlist + njoint,jointlist)          = dqddotdq;
A(jointlist + njoint,jointlist + njoint) = dqddotdqdot;
A(jointlist + njoint,end-ncontrol+1:end)   = dqddotdtau;

B = eye(D,D);
B(3,1:3) = [-1 1 1];
B(6,4:6) = [-1 1 1];
A = A * B;

V       = A';
S       = A * s * (A');
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
