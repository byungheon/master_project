%% Recursive Derivative Newton-Euler Forward Dynamics Solver
% 2018 Byungheon Kim

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  G            link inertial matrices             6*6*n
%  Vdot_0       (optional) base acceleration       6*1

%% Outputs
% [Name]        [Description]                                       [Size]
%  dtaudq       derivative of joint torques wrt joint angles        n*n
%  dtaudqdot    derivative of joint torques wrt joint velocity      n*n
%  dtaudqddot   derivative of joint torques wrt joint acceleration  n*n


%% Implementation
function [dqddotdq, dqddotdqdot, dqddotdtau] = solveForwardDynamicsDerivatives_pilco(A,M,q,qdot,qddot,G,varargin)
    %% Initialization
    n           = size(q,1);         % number of joints
    Vdot_0      = zeros(6,1);       % base acceleration
    
    if     nargin == 6         % no optional inputs
    elseif nargin == 7         % optional base acceleration
        Vdot_0 = varargin{1};        
    elseif nargin == 8
        Friction = varargin{2};
    elseif nargin == 9
        Friction = varargin{2};
        Sigmoid  = varargin{3};
    else
        error('Init Error: undefined number of inputs');
    end
    
    %%
    if nargin == 8
        [~, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, Vdot_0, Friction);
        [dtaudq, dtaudqdot, dtaudqddot] = solveInverseDynamicsDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,Vdot_0,Friction);
    elseif nargin == 9
        [~, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, Vdot_0, Friction, Sigmoid);
        [dtaudq, dtaudqdot, dtaudqddot] = solveInverseDynamicsDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,Vdot_0,Friction,Sigmoid);
    else
        [~, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, Vdot_0);
        [dtaudq, dtaudqdot, dtaudqddot] = solveInverseDynamicsDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,Vdot_0);
    end

    %% Forward Dynamics First Derivative
    dqddotdtau = pinv(dtaudqddot);
    dqddotdq    = -dqddotdtau * dtaudq;
    dqddotdqdot = -dqddotdtau * dtaudqdot;
    
end