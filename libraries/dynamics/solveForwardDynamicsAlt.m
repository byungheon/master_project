%% Closed-Form Forward Dynamics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  tau          joint torques                      n*1
%  G            link inertial matrices             6*6*n
%  Vdot_0       (optional1) base acceleration      6*1

%% Outputs
% [Name]       [Description]                      [Size]
%  qddot        joint accelerations                n*1

%% Implementation
function qddot = solveForwardDynamicsAlt(A,M,q,qdot,tau,G,varargin)
    %% Initialization
    n     = size(q,1);          % number of joints
    
    V_0    = zeros(6,1);        % base velocity
    Vdot_0 = zeros(6,1);        % base acceleration
    if     nargin == 6
    elseif nargin == 7
        Vdot_0  = varargin{1};  % optional base acceleration
    end
   %%
    M_q = zeros(n,n);
    h_q = solveInverseDynamics(A,M,q,qdot,zeros(n,1),G,Vdot_0);
    for i = 1:n
        qddot = zeros(n,1);
        qddot(i) = 1;
        M_q(:,i) = solveInverseDynamics(A,M,q,zeros(n,1),qddot,G,zeros(6,1));
    end
    %% Closed-Form Dynamics
    qddot = M_q\(tau - h_q);
end