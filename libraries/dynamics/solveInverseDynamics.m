%% Recursive Newton-Euler Inverse Dynamics Solver
% 2018 Bryan Dongik Lee

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  qddot        joint accelerations                n*1
%  G            link inertial matrices             6*6*n
%  (optional1) and (optional2)
%  Vdot_0       (optional1) base acceleration      6*1
%  F_desired    (optional2) end-effector force     6*1
%  T_end        (optional2) end-effector frame     4*4

%% Outputs
% [Name]  [Description]                      [Size]
%  tau     joint torques                      n*1
%  V       (optional) spatial velocities      6*n
%  Vdot    (optional) spatial accelerations   6*n
%  F       (optional) wrenches                6*n

%% Examples
% tau = solveInverseDynamics(A,M,q,qdot,qddot,G)
% [tau, V, Vdot] = solveInverseDynamics(A,M,q,qdot,qddot,G,Vdot_0,F_desired,T_end)
% [tau, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G,Vdot_0)

%% Implementation
function [tau, varargout] = solveInverseDynamics(A,M,q,qdot,qddot,G,varargin)
    %% Initialization
    n     = size(q,1);          % number of joints
    V     = zeros(6,n);
    Vdot  = zeros(6,n);
    F     = zeros(6,n);
    tau   = zeros(n,1);
    
    V_0    = zeros(6,1);        % base velocity
    Vdot_0 = zeros(6,1);        % base acceleration
    F_tip  = zeros(6,1);        % end-effector force
    T_end  = eye(4,4);          % end-effector pose T_{i+1,i}
    
    friction_coulomb = zeros(n,1); % coulomb friction's coefficient  
    friction_viscous = zeros(n,1); % viscous friction's coefficient 
    bSigmoid = false;
    if     nargin == 6          % no optional inputs
    elseif nargin == 7          % optional input 1
        Vdot_0 = varargin{1};  
    elseif nargin == 8
        Vdot_0           = varargin{1};         % optional base acceleration
        friction_coulomb = varargin{2}(:,1);    % optional coulomb friction's coefficient  
        friction_viscous = varargin{2}(:,2);    % optional viscous friction's coefficient   
    elseif nargin == 9
        bSigmoid         = true;
        Vdot_0           = varargin{1};         % optional base acceleration
        friction_coulomb = varargin{2}(:,1);    % optional coulomb friction's coefficient  
        friction_viscous = varargin{2}(:,2);    % optional viscous friction's coefficient 
        friction_coulomb_sigmoid = varargin{3};  % optional coulomb friction(sigmoid)'s coefficient 
    else
        error('Init Error: undefined number of inputs');
    end
    
    T    = zeros(4,4,n); % T_{i,i-1}
    Ad_T = zeros(6,6,n); % Ad_T_{i,i-1}
    
    %% Forward Recursion
    for i = 1:n
        T(:,:,i)    = exp_se3(-A(:,i)*q(i))*M(:,:,i);
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        if i == 1
            V(:,i)    = Ad_T(:,:,i)*V_0    + A(:,i)*qdot(i);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot_0 + small_ad(V(:,i))*A(:,i)*qdot(i) + A(:,i)*qddot(i);
        else
            V(:,i)    = Ad_T(:,:,i)*V(:,i-1)    + A(:,i)*qdot(i);
            Vdot(:,i) = Ad_T(:,:,i)*Vdot(:,i-1) + small_ad(V(:,i))*A(:,i)*qdot(i) + A(:,i)*qddot(i);
        end
    end

    %% Backward Recursion
    for i = n:-1:1
        if i == n
            F(:,i) = large_Ad(T_end)'*F_tip + G(:,:,i)*Vdot(:,i) - small_ad(V(:,i))'*G(:,:,i)*V(:,i);
        else
            F(:,i) = Ad_T(:,:,i+1)'*F(:,i+1) + G(:,:,i)*Vdot(:,i) - small_ad(V(:,i))'*G(:,:,i)*V(:,i);
        end
        tau(i) = A(:,i)'*F(:,i);
    end
    
    if bSigmoid
        qdot_coulomb = 2 ./ (1+exp(-friction_coulomb_sigmoid.*qdot)) - 1;
    else
        qdot_coulomb = sign(qdot);
    end
    f_c = friction_coulomb.*qdot_coulomb;
    f_v = friction_viscous.*qdot;
    tau = tau + f_c + f_v;
    if nargout > 1
        varargout{1} = T;
        varargout{2} = V;
        varargout{3} = Vdot;
        varargout{4} = F;  
    end
end