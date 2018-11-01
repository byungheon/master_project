%% gp1d.m
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

function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp1d_kuka_planar_PIREM(gpmodel, m, s)
%% Code

% If no derivatives are required, call gp1
if nargout < 4; [M S V] = gp1_kuka_planar_dyn(gpmodel, m, s); return; end
% If there are no inducing inputs, back off to gp0d (no sparse GP required)
if numel(gpmodel.induce) == 0
  [M S V dMdm dSdm dVdm dMds dSds dVds] = gp0d_kuka_planar_PIREM(gpmodel, m, s); return;
end

D_g = length(m);
gp_list = gpmodel.jointi(end)+1:length(m);
m_gp = m(gp_list);
s_gp = s(gp_list,gp_list);

% 1) If necessary: re-compute cached variables
persistent iK iK2 beta oldX;
ridge = 1e-6;                        % jitter to make matrix better conditioned
[n, D] = size(gpmodel.inputs);    % number of examples and dimension of inputs
E = size(gpmodel.targets,2);         % number of examples and number of outputs
X = gpmodel.hyp; input = gpmodel.inputs; targets = gpmodel.targets;

[np pD pE] = size(gpmodel.induce);     % number of pseudo inputs per dimension
pinput = gpmodel.induce;                                   % all pseudo inputs

if numel(X) ~= numel(oldX) || isempty(iK) || isempty(iK2) || ... % if necessary
    sum(any(X ~= oldX)) || numel(iK2) ~=E*np^2 || numel(iK) ~= n*np*E
  oldX = X;                                        % compute K, inv(K), inv(K2)
  iK = zeros(np,n,E); iK2 = zeros(np,np,E); beta = zeros(np,E);
  
  for i = 1:E
    pinp = bsxfun(@rdivide,pinput(:,:,min(i,pE)),exp(X(1:D,i)'));
    inp = bsxfun(@rdivide,input,exp(X(1:D,i)'));
    Kmm = exp(2*X(D+1,i)-maha(pinp,pinp)/2) + ridge*eye(np);  % add small ridge
    Kmn = exp(2*X(D+1,i)-maha(pinp,inp)/2);
    L = chol(Kmm)';
    V = L\Kmn;                                             % inv(sqrt(Kmm))*Kmn
    if isfield(gpmodel,'nigp')
      G = exp(2*X(D+1,i))-sum(V.^2)+gpmodel.nigp(:,i)';
    else
      G = exp(2*X(D+1,i))-sum(V.^2);
    end
    G = sqrt(1+G/exp(2*X(D+2,i)));
    V = bsxfun(@rdivide,V,G);
    Am = chol(exp(2*X(D+2,i))*eye(np) + V*V')';
    At = L*Am;                           % chol(sig*B) [Snelson's thesis, p. 40]
    iAt = At\eye(np);
    % The following is not an inverse matrix, but we'll treat it as such: multiply
    % the targets from right and the cross-covariances left to get predictive mean.
    iK(:,:,i) = ((Am\(bsxfun(@rdivide,V,G)))'*iAt)';
    beta(:,i) = iK(:,:,i)*targets(:,i);
    iB = iAt'*iAt.*exp(2*X(D+2,i));          % inv(B), [Snelson's thesis, p. 40]
    iK2(:,:,i) = Kmm\eye(np) - iB; % covariance matrix for predictive variances
  end
end

k = zeros(np,E); M = zeros(E,1); V = zeros(D,E); S = zeros(E);       % allocate
dMds = zeros(E,D,D); dSdm = zeros(E,E,D); r = zeros(1,D);
dSds = zeros(E,E,D,D); dVds = zeros(D,E,D,D); T = zeros(D);
inp = zeros(np,D,E);

% 2) Compute predicted mean and inv(s) times input-output covariance
for i = 1:E
  inp(:,:,i) = bsxfun(@minus,pinput(:,:,min(i,pE)),m_gp');   % centralize p-inputs
  
  L = diag(exp(-X(1:D,i)));
  in = inp(:,:,i)*L;
  B = L*s_gp*L+eye(D);  LiBL = L/B*L; % iR = LiBL;
  
  t = in/B;
  l = exp(-sum(in.*t,2)/2); lb = l.*beta(:,i);
  tL = t*L;
  tlb = bsxfun(@times,tL,lb);
  c = exp(2*X(D+1,i))/sqrt(det(B));
  M(i) = c*sum(lb);
  V(:,i) = tL'*lb*c;                   % inv(s) times input-output covariance
  dMds(i,:,:) = c*tL'*tlb/2 - LiBL*M(i)/2;
  for d = 1:D
    dVds(d,i,:,:) = c*bsxfun(@times,tL,tL(:,d))'*tlb/2 - LiBL*V(d,i)/2 ...
      - (V(:,i)*LiBL(d,:) + LiBL(:,d)*V(:,i)')/2;
  end
  k(:,i) = 2*X(D+1,i)-sum(in.*in,2)/2;
end
dMdm = V'; dVdm = 2*permute(dMds,[2 1 3]);                  % derivatives wrt m


iell2 = exp(-2*gpmodel.hyp(1:D,:));
inpiell2 = bsxfun(@times,inp,permute(iell2,[3,1,2])); % N-by-D-by-E

% 3) Compute predictive covariance matrix, non-central moments
for i=1:E
  ii = inpiell2(:,:,i);
  
  for j=1:i
    R = s_gp*diag(iell2(:,i)+iell2(:,j))+eye(D); t = 1/sqrt(det(R));
    ij = inpiell2(:,:,j);
    L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s_gp/2));
    
    if i==j
      iKL = iK2(:,:,i).*L; s1iKL = sum(iKL,1); s2iKL = sum(iKL,2);
      S(i,j) = t*(beta(:,i)'*L*beta(:,i) - sum(s1iKL));
      zi = ii/R;
      bibLi = L'*beta(:,i).*beta(:,i); cbLi = L'*bsxfun(@times, beta(:,i), zi);
      r = (bibLi'*zi*2 - (s2iKL' + s1iKL)*zi)*t;
      for d = 1:D
        T(d,1:d) = 2*(zi(:,1:d)'*(zi(:,d).*bibLi) + ...
          cbLi(:,1:d)'*(zi(:,d).*beta(:,i)) - zi(:,1:d)'*(zi(:,d).*s2iKL) ...
          - zi(:,1:d)'*(iKL*zi(:,d)));
        T(1:d,d) = T(d,1:d)';
      end
    else
      zi = ii/R; zj = ij/R;
      S(i,j) = beta(:,i)'*L*beta(:,j)*t; S(j,i) = S(i,j);
      
      bibLj = L*beta(:,j).*beta(:,i); bjbLi = L'*beta(:,i).*beta(:,j);
      cbLi = L'*bsxfun(@times, beta(:,i), zi);
      cbLj = L*bsxfun(@times, beta(:,j), zj);
      
      r = (bibLj'*zi+bjbLi'*zj)*t;
      for d = 1:D
        T(d,1:d) = zi(:,1:d)'*(zi(:,d).*bibLj) + ...
          cbLi(:,1:d)'*(zj(:,d).*beta(:,j)) + zj(:,1:d)'*(zj(:,d).*bjbLi) + ...
          cbLj(:,1:d)'*(zi(:,d).*beta(:,i));
        T(1:d,d) = T(d,1:d)';
      end
    end
    
    dSdm(i,j,:) = r - M(i)*dMdm(j,:)-M(j)*dMdm(i,:); dSdm(j,i,:) = dSdm(i,j,:);
    T = (t*T-S(i,j)*diag(iell2(:,i)+iell2(:,j))/R)/2;
    T = T - squeeze(M(i)*dMds(j,:,:) + M(j)*dMds(i,:,:));
    dSds(i,j,:,:) = T; dSds(j,i,:,:) = T;
  end
  
  S(i,i) = S(i,i) + exp(2*X(D+1,i));
end

% 4) Centralize moments
S = S - M*M';

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
%% converting GP variables to global
% M,S,V
invscovsx = [zeros(njoint,D);eye(D,D)];
invscovsx = sparse(invscovsx);
V = invscovsx * V;

% Gradient_g
dMdm_g = [zeros(E,njoint) dMdm];
dMds_g = zeros(E,D_g,D_g);
dMds_g(:,gp_list,gp_list) = dMds;
dSdm_g = zeros(E,E,D_g);
dSdm_g(:,:,gp_list) = dSdm;
dSds_g = zeros(E,E,D_g,D_g);
dSds_g(:,:,gp_list,gp_list) = dSds;
dVdm_g = zeros(D_g,E,D_g);
dVds_g = zeros(D_g,E,D_g,D_g);
for i = 1:D
    dVdm_g(:,:,gp_list(i)) = invscovsx * dVdm(:,:,i);
    for j = 1:D
         dVds_g(:,:,gp_list(i),gp_list(j)) = invscovsx * dVds(:,:,i,j);
    end
end
%% Robot Dynamics
n_span  = gpmodel.n_span;
D_D     = njoint*3;
dynamics_list = [jointlist njoint + jointlist [D_g-njoint+1:D_g]];
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

M(jointlist)            = M(jointlist) + (y(n_span+1,jointlist)' - q);
M(jointlist + njoint)   = M(jointlist + njoint) + (y(n_span+1,jointlist + njoint)' - qdot);

A                                        = zeros(E,D_D);
A(jointlist,jointlist + njoint)          = eye(njoint,njoint) * gpmodel.stepsize;
A(jointlist + njoint,jointlist)          = dqddotdq;
A(jointlist + njoint,jointlist + njoint) = dqddotdqdot;
A(jointlist + njoint,end-njoint+1:end)   = dqddotdtau;

%% converting Dynamics variables to global
invscovsx = zeros(D_g,D_D);
invscovsx(dynamics_list,1:end) = eye(D_D);
invscovsx = sparse(invscovsx);
V_dyn     = invscovsx * (A'); % D_g x E


%%
tmpS = V_dyn' * s * V; 
S = S + A * s(dynamics_list,dynamics_list) * (A') + tmpS + tmpS';
S = (S + S')/2;

V = V + V_dyn;

% 5) Vectorize derivatives
dMdm = dMdm_g;
dMds = reshape(dMds_g,[E D_g*D_g]);
dSds = reshape(dSds_g,[E*E D_g*D_g]); dSdm = reshape(dSdm_g,[E*E D_g]);
dVds = reshape(dVds_g,[D_g*E D_g*D_g]); dVdm = reshape(dVdm_g,[D_g*E D_g]);
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