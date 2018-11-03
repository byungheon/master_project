%% gp1d_kuka_planar_mine.m
% *Summary:* Compute joint predictions (and derivatives) for the FITC sparse
% approximation to multiple GPs with uncertain inputs.
% Predictive variances contain uncertainty about the function, but no noise.
% If gpmodel.nigp exists, individual noise contributions are added.
%
%   function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp1d(gpmodel, m, s)
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
%   dMdm       output mean by input mean                             [ E  x  D ]
%   dSdm       output covariance by input mean                       [E*E x  D ]
%   dVdm       inv(s)*input-output covariance by input mean          [D*E x  D ]
%   dMds       ouput mean by input covariance                        [ E  x D*D]
%   dSds       output covariance by input covariance                 [E*E x D*D]
%   dVds       inv(s)*input-output covariance by input covariance    [D*E x D*D]
%
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-05-24
%
%% High-Level Steps
% # If necessary, compute kernel matrix and cache it
% # Compute predicted mean and inv(s) times input-output covariance
% # Compute predictive covariance matrix, non-central moments
% # Centralize moments
% # Vectorize derivatives

function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp1d_NOGP_MINE(gpmodel, m, s)
%% Code
% If no derivatives are required, call gp1
if nargout < 4; [M S V] = gp1_NOGP(gpmodel, m, s); return; end
% If there are no inducing inputs, back off to gp0d (no sparse GP required)
% if numel(gpmodel.induce) == 0
%   [M S V dMdm dSdm dVdm dMds dSds dVds] = gp0d_NOGP_MINE(gpmodel, m, s); return;
% end
D = length(m); E = size(gpmodel.targets,2);
%% Initialization for robot dynamics
persistent dynamics OPTIONS ctrlfcn u0 par jointlist njoint ncontrol;
if isempty(dynamics)
    dynamics    = @dynamics_NOGP_learning;
    OPTIONS     = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);
    ctrlfcn     = str2func('zoh');   
    par.dt = gpmodel.stepsize; par.delay = 0; par.tau = gpmodel.stepsize;
    jointlist   = gpmodel.jointi;
    njoint      = length(jointlist);
    ncontrol    = njoint-1;
    u0          = cell(1,ncontrol);
end

%% Robot Dynamics
% n_span  = gpmodel.n_span;
q       = m(jointlist);
qdot    = m(jointlist + njoint);
tau     = [m(end-ncontrol+1:end);0];

% for j = 1:ncontrol, u0{j} = @(t)ctrlfcn(tau(j,:),t,par); end
% [~, y] = ode45(dynamics, linspace(0,gpmodel.stepsize,n_span+1), m(1:(2*njoint)), OPTIONS, u0{:});

% 
qddot   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q),stateTojoint(qdot),tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
[dqddotdq, dqddotdqdot, dqddotdtau, dqddotdqdq, dqddotdqdqdot, dqddotdqdtau, dqddotdqdotdqdot, dqddotdqdotdtau, dqddotdtaudtau] = ...
solveForwardDynamicsSecondDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q),stateTojoint(qdot),qddot,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
% for i = 2:n_span
%     q_tmp       = y(i,jointlist)';
%     qdot_tmp    = y(i,jointlist + njoint)';
%     qddot_tmp   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q_tmp),stateTojoint(qdot_tmp),tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
%     [dqddotdq_temp, dqddotdqdot_temp, dqddotdtau_temp, dqddotdqdq_temp, dqddotdqdqdot_temp, dqddotdqdtau_temp, dqddotdqdotdqdot_temp, dqddotdqdotdtau_temp, dqddotdtaudtau_temp] = ...
%         solveForwardDynamicsSecondDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,stateTojoint(q_tmp),stateTojoint(qdot_tmp),qddot_tmp,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
%     dqddotdq            = dqddotdq + dqddotdq_temp;
%     dqddotdqdot         = dqddotdqdot + dqddotdqdot_temp;
%     dqddotdtau          = dqddotdtau + dqddotdtau_temp;
%     dqddotdqdq          = dqddotdqdq + dqddotdqdq_temp;
%     dqddotdqdqdot       = dqddotdqdqdot + dqddotdqdqdot_temp;
%     dqddotdqdtau        = dqddotdqdtau + dqddotdqdtau_temp;
%     dqddotdqdotdqdot    = dqddotdqdotdqdot + dqddotdqdotdqdot_temp;
%     % no need to compute since it's zero anyway
% %     dqddotdqdotdtau     = dqddotdqdotdtau + dqddotdqdotdtau_temp; 
% %     dqddotdtaudtau      = dqddotdtaudtau + dqddotdtaudtau_temp;
% end
% dqddotdq            = dqddotdq * gpmodel.stepsize/n_span;
% dqddotdqdot         = dqddotdqdot * gpmodel.stepsize/n_span;
% dqddotdtau          = dqddotdtau(:,1:ncontrol) * gpmodel.stepsize/n_span;
dqddotdq            = dqddotdq * gpmodel.stepsize;
dqddotdqdot         = dqddotdqdot * gpmodel.stepsize;
dqddotdtau          = dqddotdtau(:,1:ncontrol) * gpmodel.stepsize;
% 
dFDdqdq         = zeros(njoint,njoint,njoint);
dFDdqdotdq      = zeros(njoint,njoint,njoint);
dFDdtaudq       = zeros(njoint,njoint,njoint);
dFDdqdqdot      = zeros(njoint,njoint,njoint);
dFDdqdotdqdot   = zeros(njoint,njoint,njoint);
dFDdtaudqdot    = zeros(njoint,njoint,njoint);
dFDdqdtau       = zeros(njoint,njoint,njoint);
dFDdqdotdtau    = zeros(njoint,njoint,njoint);
dFDdtaudtau     = zeros(njoint,njoint,njoint);
% 
% for i = 1:njoint
%     for j = 1:njoint
%         dFDdqdq(:,j,i)      = dqddotdqdq(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%         dFDdqdqdot(:,j,i)   = dqddotdqdqdot(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%         dFDdqdtau(:,j,i)    = dqddotdqdtau(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
% 
%         dFDdqdotdq(:,j,i)       = dqddotdqdqdot(:,(i-1)*njoint+j) * gpmodel.stepsize/n_span;
%         dFDdqdotdqdot(:,j,i)    = dqddotdqdotdqdot(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%         dFDdqdotdtau(:,j,i)     = dqddotdqdotdtau(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%         
%         dFDdtaudq(:,i,j)        = dqddotdqdtau(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%         dFDdtaudqdot(:,i,j)     = dqddotdqdotdtau(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%         dFDdtaudtau(:,j,i)      = dqddotdtaudtau(:,(j-1)*njoint+i) * gpmodel.stepsize/n_span;
%     end
% end

for i = 1:njoint
    for j = 1:njoint
        dFDdqdq(:,j,i)      = dqddotdqdq(:,(j-1)*njoint+i) * gpmodel.stepsize;
        dFDdqdqdot(:,j,i)   = dqddotdqdqdot(:,(j-1)*njoint+i) * gpmodel.stepsize;
        dFDdqdtau(:,j,i)    = dqddotdqdtau(:,(j-1)*njoint+i) * gpmodel.stepsize;

        dFDdqdotdq(:,j,i)       = dqddotdqdqdot(:,(i-1)*njoint+j) * gpmodel.stepsize;
        dFDdqdotdqdot(:,j,i)    = dqddotdqdotdqdot(:,(j-1)*njoint+i) * gpmodel.stepsize;
        dFDdqdotdtau(:,j,i)     = dqddotdqdotdtau(:,(j-1)*njoint+i) * gpmodel.stepsize;
        
        dFDdtaudq(:,i,j)        = dqddotdqdtau(:,(j-1)*njoint+i) * gpmodel.stepsize;
        dFDdtaudqdot(:,i,j)     = dqddotdqdotdtau(:,(j-1)*njoint+i) * gpmodel.stepsize;
        dFDdtaudtau(:,j,i)      = dqddotdtaudtau(:,(j-1)*njoint+i) * gpmodel.stepsize;
    end
end

M = zeros(E,1);
% M(jointlist)              = (y(n_span+1,jointlist)' - q);
% M(jointlist + njoint)     = (y(n_span+1,jointlist + njoint)' - qdot);
M(jointlist)                = qdot * gpmodel.stepsize;
M(jointlist + njoint)       = jointTostate(qddot) * gpmodel.stepsize;

A                                        = zeros(E,D);
A(jointlist,jointlist + njoint)          = eye(njoint,njoint) * gpmodel.stepsize;
A(jointlist + njoint,jointlist)          = dqddotdq;
A(jointlist + njoint,jointlist + njoint) = dqddotdqdot;
A(jointlist + njoint,end-ncontrol+1:end)   = dqddotdtau;

B1 = eye(D,D);
B1(3,1:3) = [-1 1 1];
B1(6,4:6) = [-1 1 1];
B2 = eye(E,E);
B2(3,1:3) = [1 -1 1];
B2(6,4:6) = [1 -1 1];

A = B2 * A * B1;

dAdm = zeros(E,D,D);
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i)                = dFDdqdq(:,:,i);
    dAdm(jointlist + njoint,jointlist + njoint,i)       = dFDdqdotdq(:,:,i);
    dAdm(jointlist + njoint,end-ncontrol+1:end,i)       = dFDdtaudq(:,1:ncontrol,i);   
end
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i+njoint)             = dFDdqdqdot(:,:,i);
    dAdm(jointlist + njoint,jointlist + njoint,i+njoint)    = dFDdqdotdqdot(:,:,i);
    dAdm(jointlist + njoint,end-ncontrol+1:end,i+njoint)      = dFDdtaudqdot(:,1:ncontrol,i);   
end
for i = 1:ncontrol
    dAdm(jointlist + njoint,jointlist,i+end-ncontrol)               = dFDdqdtau(:,:,i);
    dAdm(jointlist + njoint,jointlist + njoint,i+end-ncontrol)      = dFDdqdotdtau(:,:,i);
    dAdm(jointlist + njoint,end-ncontrol+1:end,i+end-ncontrol)      = dFDdtaudtau(:,1:ncontrol,i);   
end

dBABdA = zeros(E,D,E,D);
for i = 1:E
    for j = 1:D
        dBABdA(:,:,i,j) = B2(:,i)*B1(j,:);
    end
end
dBABdm_convert = zeros(E,D,D);
for k = 1:D
    temppp = zeros(E,D);
    for i = 1:E
        for j = 1:D
            temppp = temppp + dBABdA(:,:,i,j) * dAdm(i,j,k);
        end
    end
    dBABdm_convert(:,:,k) = temppp;
end

for i = 1:D
    temppp = zeros(E,D);
    for j = 1:D
        temppp = temppp + dBABdm_convert(:,:,j) * B1(j,i);
    end
   dAdm(:,:,i) = temppp;
end

V = A';
S = A * s * (A');

Asigma_t = A * s;

dSDdAT = zeros(E,E,D,E);
dSDdm  = zeros(E,E,D);
dSDds   = zeros(E,E,D,D);
dVDds   = zeros(D,E,D,D);

for i = 1:D
   for j = 1:E
       temp_sigma       = zeros(E,E);
       temp_sigma(:,j)  = Asigma_t(:,i);
       temp_sigma(j,:)  = temp_sigma(j,:) + Asigma_t(:,i)';
       dSDdAT(:,:,i,j)  = temp_sigma;
   end
   
   for j = 1:i
        tempmatdVD        = zeros(D,E);
        if i == j
            dSDds(:,:,i,j)    = A(:,i) * (A(:,j)');
            tempmatdVD(i,:)   = A(:,j)';
            dVDds(:,:,i,j)    = tempmatdVD;
        else
            dSDds(:,:,i,j)    = 0.5*(A(:,i) * (A(:,j)') + A(:,j) * (A(:,i)'));
            dSDds(:,:,j,i)    = dSDds(:,:,i,j);
            tempmatdVD(i,:)   = A(:,j)';
            tempmatdVD(j,:)   = tempmatdVD(j,:) + A(:,i)';
            dVDds(:,:,i,j)    = 0.5 * s \ tempmatdVD;
            dVDds(:,:,j,i)    = dVDds(:,:,i,j);
        end  
   end
end
for i = 1:D
   tempmat = zeros(E,E);
   for j = 1:D
      for k =1:E
          tempmat = tempmat + dSDdAT(:,:,j,k) * dAdm(k,j,i);
      end
   end
   dSDdm(:,:,i)     = tempmat;
end

dMdm = A;
dMds = zeros(E,D,D);
dSdm = dSDdm;
dSds = dSDds;
dVds = dVDds;
dVdm = permute(dAdm,[2,1,3]);
%%

% 5) vectorize derivatives
% dMds = reshape(dMds,[E D*D]);
% % dSds = kron(A,A); % kron(A_g,A_g) == dSDds
% dSds = reshape(dSds,[E*E,D*D]);
% dSdm = reshape(dSdm,[E*E D]);
% dVds = reshape(dVds,[D*E D*D]);
% dVdm = reshape(dVdm,[D*E D]);

% 6) symmetrize variables dependent to variance
% X=reshape(1:D*D,[D D]); XT=X'; dSds=(dSds+dSds(:,XT(:)))/2; 
% dMds=(dMds+dMds(:,XT(:)))/2;
% dVds=(dVds+dVds(:,XT(:)))/2;
% X=reshape(1:E*E,[E E]); XT=X'; dSds=(dSds+dSds(XT(:),:))/2;
% dSdm=(dSdm+dSdm(XT(:),:))/2; 

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