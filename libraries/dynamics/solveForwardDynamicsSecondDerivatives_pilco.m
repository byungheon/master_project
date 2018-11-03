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
function [dqddotdq, dqddotdqdot, dqddotdtau, dqddotdqdq, dqddotdqdqdot, dqddotdqdtau, dqddotdqdotdqdot, dqddotdqdotdtau, dqddotdtaudtau] = solveForwardDynamicsSecondDerivatives_pilco(A,M,q,qdot,qddot,G,varargin)
    %% Initialization
    n           = size(q,1);         % number of joints
    Vdot_0      = zeros(6,1);       % base acceleration

    if     nargin == 6         % no optional inputs
    elseif nargin == 7         % optional base acceleration
        Vdot_0 = varargin{1};        
    elseif nargin == 8
        Vdot_0 = varargin{1}; 
        Friction = varargin{2};
    elseif nargin == 9
        Vdot_0 = varargin{1}; 
        Friction = varargin{2};
        Sigmoid  = varargin{3};
    else
        error('Init Error: undefined number of inputs');
    end
    
    %%
    if nargin == 8
        [~, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, Vdot_0, Friction);
        [dtaudq, dtaudqdot, dtaudqddot, dtaudqdq, dtaudqdqdot, dtaudqdqddot, dtaudqdotdqdot, dtaudqdotdqddot, dtaudqddotdqddot] = solveInverseDynamicsSecondDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,Vdot_0,Friction);
    elseif nargin == 9
        [~, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, Vdot_0, Friction, Sigmoid);
        [dtaudq, dtaudqdot, dtaudqddot, dtaudqdq, dtaudqdqdot, dtaudqdqddot, dtaudqdotdqdot, dtaudqdotdqddot, dtaudqddotdqddot] = solveInverseDynamicsSecondDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,Vdot_0,Friction,Sigmoid);
    else
        [~, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, Vdot_0);
        [dtaudq, dtaudqdot, dtaudqddot, dtaudqdq, dtaudqdqdot, dtaudqdqddot, dtaudqdotdqdot, dtaudqdotdqddot, dtaudqddotdqddot] = solveInverseDynamicsSecondDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,Vdot_0);
    end
    
    
    %% Forward Dynamics First Derivative
    dqddotdtau = pinv(dtaudqddot);
    dqddotdq    = -dtaudqddot \ dtaudq;
    dqddotdqdot = -dtaudqddot \ dtaudqdot;
    %% Forward Dynamics Second Derivative
    dqddotdqdq          = zeros(n, n*n);
    dqddotdqdqdot       = zeros(n, n*n);
    dqddotdqdtau        = zeros(n, n*n);
    dqddotdqdotdqdot    = zeros(n, n*n);
    dqddotdqdotdtau     = zeros(n, n*n);
    dqddotdtaudtau      = zeros(n, n*n);
    
    for i = 1:n
       for j = 1:n
           dqddotdqdq(:,(i-1) * n + j)          = dtaudqdq(:,(i-1) * n + j);
           dqddotdqdqdot(:,(i-1) * n + j)       = dtaudqdqdot(:,(i-1) * n + j);
           dqddotdqdtau(:,(i-1) * n + j)        = zeros(n,1);
           dqddotdqdotdqdot(:,(i-1) * n + j)    = dtaudqdotdqdot(:,(i-1) * n + j);
           for k = 1:n
              if j == i
                  dqddotdqdq(:,(i-1) * n + j)       = dqddotdqdq(:,(i-1) * n + j) + 2 * dtaudqdqddot(:,(i-1) * n + k) * dqddotdq(k, j);
              else
                  dqddotdqdq(:,(i-1) * n + j)       = dqddotdqdq(:,(i-1) * n + j) + dtaudqdqddot(:,(i-1) * n + k) * dqddotdq(k, j)  + dtaudqdqddot(:,(j-1) * n + k) * dqddotdq(k, i);
              end

              dqddotdqdqdot(:,(i-1) * n + j)    = dqddotdqdqdot(:,(i-1) * n + j) + dtaudqdqddot(:,(i-1) * n + k) * dqddotdqdot(k, j);
              dqddotdqdtau(:,(i-1) * n + j)     = dqddotdqdtau(:,(i-1) * n + j) + dtaudqdqddot(:,(i-1) * n + k) * dqddotdtau(k, j); 
           end
           dqddotdqdq(:,(i-1) * n + j)          = -dtaudqddot \ dqddotdqdq(:,(i-1) * n + j);
           dqddotdqdqdot(:,(i-1) * n + j)       = -dtaudqddot \ dqddotdqdqdot(:,(i-1) * n + j);
           dqddotdqdtau(:,(i-1) * n + j)        = -dtaudqddot \ dqddotdqdtau(:,(i-1) * n + j);
           dqddotdqdotdqdot(:,(i-1) * n + j)    = -dtaudqddot \ dqddotdqdotdqdot(:,(i-1) * n + j);
       end
    end
end