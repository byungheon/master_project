%% gp0d_kuka_mine.m
% *Summary:* Compute joint predictions and derivatives for multiple GPs
% with uncertain inputs. Predictive variances contain uncertainty about the 
% function, but no noise.
% If gpmodel.nigp exists, individial noise contributions are added.
% 
%
%   function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp0d_kuka_mine(gpmodel, m, s)
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

function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp0d_kuka_planar_mine(gpmodel, m, s)
%% Code 
% If no derivatives required, call gp0
if nargout < 4; [M S V] = gp0_kuka_dyn(gpmodel, m, s); return; end

%% GP
persistent K iK beta oldX oldn;
[n, D] = size(gpmodel.inputs);    % number of examples and dimension of inputs
E = size(gpmodel.targets,2);                               % number of outputs
X = gpmodel.hyp;                              % short hand for hyperparameters

% 1) if necessary: re-compute cached variables
if numel(X) ~= numel(oldX) || isempty(iK) || sum(any(X ~= oldX)) || n ~= oldn
  oldX = X; oldn = n;
  iK = zeros(n,n,E);  K = iK; beta = zeros(n,E);
  
  for i=1:E                                              % compute K and inv(K)
    inp = bsxfun(@rdivide,gpmodel.inputs,exp(X(1:D,i)'));
    K(:,:,i) = exp(2*X(D+1,i)-maha(inp,inp)/2);
    if isfield(gpmodel,'nigp')
      L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n) + diag(gpmodel.nigp(:,i)))';
    else
      L = chol(K(:,:,i) + exp(2*X(D+2,i))*eye(n))';
    end
    iK(:,:,i) = L'\(L\eye(n));
    beta(:,i) = L'\(L\gpmodel.targets(:,i));
  end
end

k = zeros(n,E); M = zeros(E,1); V = zeros(D,E); S = zeros(E);      % initialize
dMds = zeros(E,D,D); dSdm = zeros(E,E,D);
dSds = zeros(E,E,D,D); dVds = zeros(D,E,D,D); T = zeros(D);

inp = bsxfun(@minus,gpmodel.inputs,m');                     % centralize inputs

% 2) compute predicted mean and inv(s) times input-output covariance
for i=1:E
  iL = diag(exp(-X(1:D,i))); % inverse length scales
  in = inp*iL;
  B = iL*s*iL+eye(D); LiBL = iL/B*iL;
  t = in/B;
  l = exp(-sum(in.*t,2)/2); lb = l.*beta(:,i);
  tL = t*iL;
  tlb = bsxfun(@times,tL,lb);
  c = exp(2*X(D+1,i))/sqrt(det(B));
  M(i) = c*sum(lb);
  V(:,i) = tL'*lb*c;                     % inv(s) times input-output covariance
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

% 3) compute predictive covariance matrix, non-central moments
for i=1:E
  ii = inpiell2(:,:,i);
  
  for j=1:i
    R = s*diag(iell2(:,i)+iell2(:,j))+eye(D); 
    t = 1/sqrt(det(R));
    ij = inpiell2(:,:,j);
    L = exp(bsxfun(@plus,k(:,i),k(:,j)')+maha(ii,-ij,R\s/2));
    if i==j
      iKL = iK(:,:,i).*L; 
      s1iKL = sum(iKL,1); 
      s2iKL = sum(iKL,2);
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
      S(i,j) = beta(:,i)'*L*beta(:,j)*t; 
      S(j,i) = S(i,j);
      
      bibLj = L*beta(:,j).*beta(:,i); 
      bjbLi = L'*beta(:,i).*beta(:,j);
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
    
    dSdm(i,j,:) = r - M(i)*dMdm(j,:)-M(j)*dMdm(i,:); 
    dSdm(j,i,:) = dSdm(i,j,:);
    T = (t*T-S(i,j)*diag(iell2(:,i)+iell2(:,j))/R)/2;
    T = T - reshape(M(i)*dMds(j,:,:) + M(j)*dMds(i,:,:),D,D);
    dSds(i,j,:,:) = T; 
    dSds(j,i,:,:) = T;
  end
  
  S(i,i) = S(i,i) + exp(2*X(D+1,i));
end
% 4) centralize moments
S = S - M*M';
%% Robot Dynamics
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

q       = m(jointlist);
qdot    = m(jointlist + njoint);
tau     = m(end-njoint+1:end);

for j = 1:njoint, u0{j} = @(t)ctrlfcn(tau(j,:),t,par); end
[~, y] = ode45(dynamics, [0 gpmodel.stepsize/2 gpmodel.stepsize], m(1:(2*njoint)), OPTIONS, u0{:});

% qddot   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,q,qdot,tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
qddot = (y(3,jointlist + njoint)' - qdot)/gpmodel.stepsize;
[dqddotdq, dqddotdqdot, dqddotdtau, dqddotdqdq, dqddotdqdqdot, dqddotdqdtau, dqddotdqdotdqdot, dqddotdqdotdtau, dqddotdtaudtau] = ...
solveForwardDynamicsSecondDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,q,qdot,qddot,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
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

for i = 1:njoint
    for j = 1:njoint
        dFDdqdq(:,j,i)      = dqddotdqdq(:,(j-1)*njoint+i);
        dFDdqdqdot(:,j,i)   = dqddotdqdqdot(:,(j-1)*njoint+i);
        dFDdqdtau(:,j,i)    = dqddotdqdtau(:,(j-1)*njoint+i);

        dFDdqdotdq(:,j,i)       = dqddotdqdqdot(:,(i-1)*njoint+j);
        dFDdqdotdqdot(:,j,i)    = dqddotdqdotdqdot(:,(j-1)*njoint+i);
        dFDdqdotdtau(:,j,i)     = dqddotdqdotdtau(:,(j-1)*njoint+i);
        
        dFDdtaudq(:,i,j)        = dqddotdqdtau(:,(j-1)*njoint+i);
        dFDdtaudqdot(:,i,j)     = dqddotdqdotdtau(:,(j-1)*njoint+i);
        dFDdtaudtau(:,j,i)      = dqddotdtaudtau(:,(j-1)*njoint+i);
    end
end



M(jointlist)            = M(jointlist) + (y(3,jointlist)' - q);
M(jointlist + njoint)   = M(jointlist + njoint) + (y(3,jointlist + njoint)' - qdot);

A                                        = zeros(E,D);
A(jointlist,jointlist + njoint)          = eye(njoint,njoint) * gpmodel.stepsize;
A(jointlist + njoint,jointlist)          = dqddotdq * gpmodel.stepsize;
A(jointlist + njoint,jointlist + njoint) = dqddotdqdot * gpmodel.stepsize;
A(jointlist + njoint,end-njoint+1:end)   = dqddotdtau * gpmodel.stepsize;

V = V + A';

Asigma_t = A * s;

S = S + Asigma_t * V; 

cov = s * V;

dAdm = zeros(E,D,D);
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i)                = dFDdqdq(:,:,i) * gpmodel.stepsize;
    dAdm(jointlist + njoint,jointlist + njoint,i)       = dFDdqdotdq(:,:,i) * gpmodel.stepsize;
    dAdm(jointlist + njoint,end-njoint+1:end,i)         = dFDdtaudq(:,:,i) * gpmodel.stepsize;   
end
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i+njoint)             = dFDdqdqdot(:,:,i) * gpmodel.stepsize;
    dAdm(jointlist + njoint,jointlist + njoint,i+njoint)    = dFDdqdotdqdot(:,:,i) * gpmodel.stepsize;
    dAdm(jointlist + njoint,end-njoint+1:end,i+njoint)      = dFDdtaudqdot(:,:,i) * gpmodel.stepsize;   
end
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i+end-njoint)             = dFDdqdtau(:,:,i) * gpmodel.stepsize;
    dAdm(jointlist + njoint,jointlist + njoint,i+end-njoint)    = dFDdqdotdtau(:,:,i) * gpmodel.stepsize;
    dAdm(jointlist + njoint,end-njoint+1:end,i+end-njoint)      = dFDdtaudtau(:,:,i) * gpmodel.stepsize;   
end
dSDdAT = zeros(E,E,D,E);
dSDdm  = zeros(E,E,D);
for i = 1:D
   for j = 1:E
       temp_sigma       = zeros(E,E);
       temp_sigma(:,j)  = Asigma_t(:,i);
       temp_sigma(j,:)  = temp_sigma(j,:) + Asigma_t(:,i)';
       dSDdAT(:,:,i,j)  = temp_sigma;
   end
end
dAcovdm = zeros(E,E,D);
Adcovds = zeros(E,E,D,D);
dSDds   = zeros(E,E,D,D);
dVDds   = zeros(D,E,D,D);
for i = 1:D
   tempmat = zeros(E,E);
   for j = 1:D
      tempmatdVD        = zeros(D,E);
      tempmatdVD(i,:)   = A(:,j)';
      dSDds(:,:,i,j)    = A(:,i) * (A(:,j)');
      Adcovds(:,:,i,j)  = Asigma_t * dVds(:,:,i,j);  
      dVDds(:,:,i,j)    = s \ tempmatdVD;
      for k =1:E
          tempmat = tempmat + dSDdAT(:,:,j,k) * dAdm(k,j,i);
      end
   end
   dSDdm(:,:,i)     = tempmat;
   dAcovdm(:,:,i)   = dAdm(:,:,i) * cov + Asigma_t * dVdm(:,:,i);
end

dMdm = dMdm + A;
dSdm = dSdm + dSDdm + dAcovdm;
dSds = dSds + dSDds + Adcovds;
dVdm = dVdm + permute(dAdm,[2,1,3]);
dVds = dVds + dVDds;
%%
% 5) vectorize derivatives
dMds = reshape(dMds,[E D*D]);
dSds = reshape(dSds,[E*E D*D]); dSdm = reshape(dSdm,[E*E D]);
dVds = reshape(dVds,[D*E D*D]); dVdm = reshape(dVdm,[D*E D]);
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