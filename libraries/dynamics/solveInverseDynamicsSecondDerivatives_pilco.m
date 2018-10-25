%% Recursive Second Derivative Newton-Euler Inverse Dynamics Solver
% 2018 Byungheon Kim

%% Inputs
% [Name]       [Description]                      [Size]
%  A            i-th body screws from i-th frame   6*n    n = num of joints
%  M            initial relative frames M_{i,i-1}  4*4*n
%  q            joint angles                       n*1
%  qdot         joint vleocities                   n*1
%  G            link inertial matrices             6*6*n
%  T            relative link frames               4*4*n
%  V            spatial velocities                 6*n
%  Vdot         spatial accelerations              6*n
%  F            wrenches                           6*n
%  Vdot_0       (optional) base acceleration       6*1

%% Outputs
% [Name]            [Description]                                       [Size]
%  dtaudq           derivative of joint torques wrt joint angles        n*n
%  dtaudqdot        derivative of joint torques wrt joint velocity      n*n
%  dtaudqddot       derivative of joint torques wrt joint acceleration  n*n
%  dtaudqdq         second derivative of joint torques wrt q, q         n*n^2
%  dtaudqdqdot      second derivative of joint torques wrt q, qdot      n*n^2
%  dtaudqdqddot     second derivative of joint torques wrt q, qddot     n*n^2
%  dtaudqdotdqdot   second derivative of joint torques wrt qdot, qdot   n*n^2
%  dtaudqdotdqddot  second derivative of joint torques wrt qdot, qddot  n*n^2
%  dtaudqddotdqddot second derivative of joint torques wrt qddot, qddot n*n^2

%% Implementation
function [dtaudq, dtaudqdot, dtaudqddot, dtaudqdq, dtaudqdqdot, dtaudqdqddot, dtaudqdotdqdot, dtaudqdotdqddot, dtaudqddotdqddot] = solveInverseDynamicsSecondDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,varargin)
    %% Initialization
    n           = size(q,1);         % number of joints
   
    dVdq        = zeros(6,n,n);      % dV x dq x i(link index)
    dVdqdot     = zeros(6,n,n);              
    dVdotdq     = zeros(6,n,n);
    dVdotdqdot  = zeros(6,n,n);
    dVdotdqddot = zeros(6,n,n);
    
    dVdqdq          = zeros(6,n*n,n);
    dVdqdqdot   	= zeros(6,n*n,n);
    dVdotdqdq       = zeros(6,n*n,n);
    dVdotdqdqdot    = zeros(6,n*n,n);
    dVdotdqdqddot   = zeros(6,n*n,n);
    dVdotdqdotdqdot = zeros(6,n*n,n);
    
%     dFdqtmp         = zeros(6,n);
%     dFdqdottmp      = zeros(6,n);
%     dFdqddottmp     = zeros(6,n);
    
    dFdq            = zeros(6,n,n);
    dFdqdot         = zeros(6,n,n);
    dFdqddot        = zeros(6,n,n);
    
    dFdqdq          = zeros(6,n*n,n);
    dFdqdqdot   	= zeros(6,n*n,n);
    dFdqdqddot      = zeros(6,n*n,n);
    dFdqdotdqdot    = zeros(6,n*n,n);
    
    dtaudq          = zeros(n,n);
    dtaudqdot       = zeros(n,n);
    dtaudqddot      = zeros(n,n);
    
    dtaudqdq          = zeros(n,n*n);
    dtaudqdqdot       = zeros(n,n*n);
    dtaudqdqddot      = zeros(n,n*n);
    dtaudqdotdqdot    = zeros(n,n*n);
    dtaudqdotdqddot   = zeros(n,n*n);
    dtaudqddotdqddot  = zeros(n,n*n);
    
    V_0     = zeros(6,1);       % base velocity
    Vdot_0  = zeros(6,1);       % base acceleration
    % derivative of base velocity & acceleration
%     dVdq_0          = zeros(6,n);
%     dVdqdot_0       = zeros(6,n);
%     dVdotdq_0       = zeros(6,n);
%     dVdotdqdot_0    = zeros(6,n);
%     dVdotdqddot_0   = zeros(6,n);
    bSigmoid = false;
    friction_coulomb = zeros(n,1); % coulomb friction's coefficient  
    friction_viscous = zeros(n,1); % viscous friction's coefficient 
    if     nargin == 9         % no optional inputs
    elseif nargin == 10         % optional base acceleration
        Vdot_0 = varargin{1};        
    elseif nargin == 11
        Vdot_0           = varargin{1};         % optional base acceleration
        friction_coulomb = varargin{2}(:,1);    % optional coulomb friction's coefficient  
        friction_viscous = varargin{2}(:,2);    % optional viscous friction's coefficient   
    elseif nargin == 12
        bSigmoid         = true;
        Vdot_0           = varargin{1};         % optional base acceleration
        friction_coulomb = varargin{2}(:,1);    % optional coulomb friction's coefficient  
        friction_viscous = varargin{2}(:,2);    % optional viscous friction's coefficient 
        friction_coulomb_sigmoid = varargin{3};  % optional coulomb friction(sigmoid)'s coefficient 
    else
        error('Init Error: undefined number of inputs');
    end
    
    Ad_T = zeros(6,6,n); % Ad_T_{i,i-1}
    
    %% Forward Recursion
    for i = 1:n
        Ad_T(:,:,i) = large_Ad(T(:,:,i));
        ad_A = small_ad(A(:,i));
        if i == 1
%             dVdq(:,:,i)     = Ad_T(:,:,i) * dVdq_0;
            dVdq(:,i,i)     = dVdq(:,i,i) - ad_A * Ad_T(:,:,i) * V_0;
%             dVdqdot(:,:,i)  = Ad_T(:,:,i) * dVdqdot_0;
            dVdqdot(:,i,i)  = dVdqdot(:,i,i) + A(:,i);
         
%             dVdotdq(:,:,i)  = Ad_T(:,:,i) * dVdotdq_0 - ad_A * dVdq(:,:,i) * qdot(i);
            dVdotdq(:,:,i)  = - ad_A * dVdq(:,:,i) * qdot(i);
            dVdotdq(:,i,i)  = dVdotdq(:,i,i) - ad_A * Ad_T(:,:,i) * Vdot_0;
%             dVdotdqdot(:,:,i) = Ad_T(:,:,i) * dVdotdqdot_0 - ad_A * dVdqdot(:,:,i) * qdot(i);
            dVdotdqdot(:,:,i) = - ad_A * dVdqdot(:,:,i) * qdot(i);
            dVdotdqdot(:,i,i) = dVdotdqdot(:,i,i) - ad_A * V(:,i);
%             dVdotdqddot(:,:,i) = Ad_T(:,:,i) * dVdotdqddot_0;
            dVdotdqddot(:,i,i) = dVdotdqddot(:,i,i) + A(:,i);

            for j=1:n
               for k = 1:n
                  dVdotdqdq(:,(j-1)*n+k,i) = dVdotdqdq(:,(j-1)*n+k,i) + small_ad(dVdqdq(:,(j-1)*n+k,i)) * A(:,i) * qdot(i);
                  if k == i
                     dVdotdqdqdot(:,(j-1)*n+k,i) = dVdotdqdqdot(:,(j-1)*n+k,i) - ad_A * dVdq(:,j,i);
                  end
                  dVdotdqdqdot(:,(j-1)*n+k,i) = dVdotdqdqdot(:,(j-1)*n+k,i) + small_ad(dVdqdqdot(:,(j-1)*n+k,i)) * A(:,i) * qdot(i);
                  if j == i                     
                     dVdotdqdotdqdot(:,(j-1)*n+k,i) = dVdotdqdotdqdot(:,(j-1)*n+k,i) - ad_A * dVdqdot(:,k,i);
                  end
                  if k == i
                     dVdotdqdotdqdot(:,(j-1)*n+k,i) = dVdotdqdotdqdot(:,(j-1)*n+k,i) - ad_A * dVdqdot(:,j,i);
                  end
               end
            end
            dVdotdqdq(:,(i-1)+i,i) = dVdotdqdq(:,(i-1)+i,i) + ad_A * ad_A * Ad_T(:,:,i) * Vdot_0;
        else
            dVdq(:,:,i)     = Ad_T(:,:,i) * dVdq(:,:,i-1);
            dVdq(:,i,i)     = dVdq(:,i,i) - ad_A * Ad_T(:,:,i) * V(:,i-1);
            dVdqdot(:,:,i)  = Ad_T(:,:,i) * dVdqdot(:,:,i-1);
            dVdqdot(:,i,i)  = dVdqdot(:,i,i) + A(:,i);
            
            for j=1:n*n
               dVdqdq(:,j,i)    = Ad_T(:,:,i) * dVdqdq(:,j,i-1);
               dVdqdqdot(:,j,i) = Ad_T(:,:,i) * dVdqdqdot(:,j,i-1);
            end
            
            for j = 1: n
               dVdqdq(:,(i-1)*n+j,i)     = dVdqdq(:,(i-1)*n+j,i) - ad_A * Ad_T(:,:,i) * dVdq(:,j,i-1);
               dVdqdq(:,(j-1) * n + i,i) = dVdqdq(:,(j-1) * n + i,i) - ad_A * Ad_T(:,:,i) * dVdq(:,j,i-1);
               dVdqdqdot(:,(i-1)*n+j,i)  = dVdqdqdot(:,(i-1)*n+j,i) - ad_A * Ad_T(:,:,i) * dVdqdot(:,j,i-1);
            end
            dVdqdq(:,(i-1)*n+i,i) = dVdqdq(:,(i-1)*n+i,i) + ad_A * ad_A * Ad_T(:,:,i) * V(:,i-1);
            
            dVdotdq(:,:,i)  = Ad_T(:,:,i) * dVdotdq(:,:,i-1) - ad_A * dVdq(:,:,i) * qdot(i);
            dVdotdq(:,i,i)  = dVdotdq(:,i,i) - ad_A * Ad_T(:,:,i) * Vdot(:,i-1);
            dVdotdqdot(:,:,i) = Ad_T(:,:,i) * dVdotdqdot(:,:,i-1) - ad_A * dVdqdot(:,:,i) * qdot(i);
            dVdotdqdot(:,i,i) = dVdotdqdot(:,i,i) - ad_A * V(:,i);
            dVdotdqddot(:,:,i) = Ad_T(:,:,i) * dVdotdqddot(:,:,i-1);
            dVdotdqddot(:,i,i) = dVdotdqddot(:,i,i) + A(:,i);
            
            for j=1:n*n
               dVdotdqdq(:,j,i)         = Ad_T(:,:,i) * dVdotdqdq(:,j,i-1);
               dVdotdqdqdot(:,j,i)      = Ad_T(:,:,i) * dVdotdqdqdot(:,j,i-1);
               dVdotdqdqddot(:,j,i)     = Ad_T(:,:,i) * dVdotdqdqddot(:,j,i-1);
               dVdotdqdotdqdot(:,j,i)   = Ad_T(:,:,i) * dVdotdqdotdqdot(:,j,i-1);
            end
            
            for j=1:n
               for k = 1:n
                  if j == i
                     dVdotdqdq(:,(j-1)*n+k,i) = dVdotdqdq(:,(j-1)*n+k,i) - ad_A * Ad_T(:,:,i) * dVdotdq(:,k,i-1);
                  end
                  if k == i
                     dVdotdqdq(:,(j-1)*n+k,i) = dVdotdqdq(:,(j-1)*n+k,i) - ad_A * Ad_T(:,:,i) * dVdotdq(:,j,i-1);
                  end
                  dVdotdqdq(:,(j-1)*n+k,i) = dVdotdqdq(:,(j-1)*n+k,i) + small_ad(dVdqdq(:,(j-1)*n+k,i)) * A(:,i) * qdot(i);
                  
                  if j == i
                     dVdotdqdqdot(:,(j-1)*n+k,i) = dVdotdqdqdot(:,(j-1)*n+k,i) - ad_A * Ad_T(:,:,i) * dVdotdqdot(:,k,i-1);
                  end
                  if k == i
                     dVdotdqdqdot(:,(j-1)*n+k,i) = dVdotdqdqdot(:,(j-1)*n+k,i) - ad_A * dVdq(:,j,i);
                  end
                  dVdotdqdqdot(:,(j-1)*n+k,i) = dVdotdqdqdot(:,(j-1)*n+k,i) + small_ad(dVdqdqdot(:,(j-1)*n+k,i)) * A(:,i) * qdot(i);
                  if j == i
                     dVdotdqdqddot(:,(j-1)*n+k,i) = dVdotdqdqddot(:,(j-1)*n+k,i) - ad_A * Ad_T(:,:,i) * dVdotdqddot(:,k,i-1);
                     
                     dVdotdqdotdqdot(:,(j-1)*n+k,i) = dVdotdqdotdqdot(:,(j-1)*n+k,i) - ad_A * dVdqdot(:,k,i);
                  end
                  if k == i
                     dVdotdqdotdqdot(:,(j-1)*n+k,i) = dVdotdqdotdqdot(:,(j-1)*n+k,i) - ad_A * dVdqdot(:,j,i);
                  end
               end
            end
            dVdotdqdq(:,(i-1)*n+i,i) = dVdotdqdq(:,(i-1)*n+i,i) + ad_A * ad_A * Ad_T(:,:,i) * Vdot(:,i-1);
        end
    end

    %% Backward Recursion
    for i = n:-1:1
        ad_A             = small_ad(A(:,i));
        dFdqtmp          = dFdq(:,:,i) + G(:,:,i) * dVdotdq(:,:,i) - small_ad(V(:,i))' * G(:,:,i) * dVdq(:,:,i);
        dFdqdottmp       = dFdqdot(:,:,i) + G(:,:,i) * dVdotdqdot(:,:,i) - small_ad(V(:,i))' * G(:,:,i) * dVdqdot(:,:,i);
        dFdqddottmp      = dFdqddot(:,:,i) + G(:,:,i) * dVdotdqddot(:,:,i);
        
        tmp = G(:,:,i) * V(:,i);
        for j = 1:n
           dFdqtmp(:,j)     = dFdqtmp(:,j) - small_ad(dVdq(:,j,i))' * tmp;
           dFdqdottmp(:,j)  = dFdqdottmp(:,j) - small_ad(dVdqdot(:,j,i))' * tmp;
        end
        
        dtaudq(i,:)     = A(:,i)' * dFdqtmp;
        dtaudqdot(i,:)  = A(:,i)' * dFdqdottmp;
        dtaudqddot(i,:) = A(:,i)' * dFdqddottmp;
        
        dFdqdqtmp       = zeros(6, n*n);
        dFdqdqdottmp    = zeros(6, n*n);
        dFdqdqddottmp   = zeros(6, n*n);
        dFdqdotdqdottmp = zeros(6, n*n);
        
        for j = 1:n
           for k = 1:n                                                     
              dFdqdqtmp(:,(j-1)*n+k)        = dFdqdq(:,(j-1)*n+k,i) + G(:,:,i) * dVdotdqdq(:,(j-1)*n+k,i) - small_ad(V(:,i))' * G(:,:,i) * dVdqdq(:,(j-1)*n+k,i) ...
                                            - small_ad(dVdq(:,j,i))' * G(:,:,i) * dVdq(:,k,i) - small_ad(dVdq(:,k,i))' * G(:,:,i) * dVdq(:,j,i) ...  
                                            - small_ad(dVdqdq(:,(j-1)*n+k,i))' * G(:,:,i) * V(:,i);
              dFdqdqdottmp(:,(j-1)*n+k)     = dFdqdqdot(:,(j-1)*n+k,i) + G(:,:,i) * dVdotdqdqdot(:,(j-1)*n+k,i) - small_ad(V(:,i))' * G(:,:,i) * dVdqdqdot(:,(j-1)*n+k,i) ...
                                            - small_ad(dVdq(:,j,i))' * G(:,:,i) * dVdqdot(:,k,i) - small_ad(dVdqdot(:,k,i))' * G(:,:,i) * dVdq(:,j,i) ...  
                                            - small_ad(dVdqdqdot(:,(j-1)*n+k,i))' * G(:,:,i) * V(:,i);
              dFdqdqddottmp(:,(j-1)*n+k)    = dFdqdqddot(:,(j-1)*n+k,i) + G(:,:,i) * dVdotdqdqddot(:,(j-1)*n+k,i);
              dFdqdotdqdottmp(:,(j-1)*n+k)  = dFdqdotdqdot(:,(j-1)*n+k,i) + G(:,:,i) * dVdotdqdotdqdot(:,(j-1)*n+k,i) ...
                                            - small_ad(dVdqdot(:,j,i))' * G(:,:,i) * dVdqdot(:,k,i) - small_ad(dVdqdot(:,k,i))' * G(:,:,i) * dVdqdot(:,j,i);
           end
        end
        if i~=1
           dFdq(:,:,i-1)    = dFdq(:,:,i-1) + Ad_T(:,:,i)' * dFdqtmp;
           dFdq(:,i,i-1)    = dFdq(:,i,i-1) + Ad_T(:,:,i)' * small_ad(-A(:,i))' * F(:,i);
           dFdqdot(:,:,i-1) = dFdqdot(:,:,i-1) + Ad_T(:,:,i)' * dFdqdottmp;
           dFdqddot(:,:,i-1)= dFdqddot(:,:,i-1) + Ad_T(:,:,i)' * dFdqddottmp;
        end
        
        for j=1:n*n
           dtaudqdq(i,j)        = A(:,i)'*dFdqdqtmp(:,j);
           dtaudqdqdot(i,j)     = A(:,i)'*dFdqdqdottmp(:,j);
           dtaudqdqddot(i,j)    = A(:,i)'*dFdqdqddottmp(:,j);
           dtaudqdotdqdot(i,j)  = A(:,i)'*dFdqdotdqdottmp(:,j);
           if i~=1
               dFdqdq(:,j,i-1)       = dFdqdq(:,j,i-1) + Ad_T(:,:,i)' * dFdqdqtmp(:,j);
               dFdqdqdot(:,j,i-1)    = dFdqdqdot(:,j,i-1) + Ad_T(:,:,i)' * dFdqdqdottmp(:,j);
               dFdqdqddot(:,j,i-1)   = dFdqdqddot(:,j,i-1) + Ad_T(:,:,i)' * dFdqdqddottmp(:,j);
               dFdqdotdqdot(:,j,i-1) = dFdqdotdqdot(:,j,i-1) + Ad_T(:,:,i)' * dFdqdotdqdottmp(:,j);
           end
        end
        if i~=1
            for j=1:n
                dFdqdq(:,(i-1)*n+j,i-1)     = dFdqdq(:,(i-1)*n+j,i-1) + Ad_T(:,:,i)' * small_ad(-A(:,i))' * dFdqtmp(:,j);
                dFdqdq(:,(j-1)*n+i,i-1)     = dFdqdq(:,(j-1)*n+i,i-1) + Ad_T(:,:,i)' * small_ad(-A(:,i))' * dFdqtmp(:,j);
                dFdqdqdot(:,(i-1)*n+j,i-1)  = dFdqdqdot(:,(i-1)*n+j,i-1) + Ad_T(:,:,i)' * small_ad(-A(:,i))' * dFdqdottmp(:,j);
                dFdqdqddot(:,(i-1)*n+j,i-1) = dFdqdqddot(:,(i-1)*n+j,i-1) + Ad_T(:,:,i)' * small_ad(-A(:,i))' * dFdqddottmp(:,j);
            end
            dFdqdq(:,(i-1)*n+i,i-1) = dFdqdq(:,(i-1)*n+i,i-1) + Ad_T(:,:,i)' * ad_A' * ad_A' * F(:,i);
        end
    end
    if bSigmoid 
        exp_sig = exp(-friction_coulomb_sigmoid.*qdot);
        exp_sig_2 = (1+exp_sig).^2;
        dtaudq_friction = diag(2*friction_coulomb.*friction_coulomb_sigmoid.*exp_sig./exp_sig_2 + friction_viscous);
    else
        dtaudq_friction = diag(friction_viscous);
    end
    dtaudqdot = dtaudqdot + dtaudq_friction;
end