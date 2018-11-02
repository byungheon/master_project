
function [M, S, V, dMdm, dSdm, dVdm, dMds, dSds, dVds] = gp1d_odetest(gpmodel, m, s)
%% Code
% if ~isfield(gpmodel,'induce') || numel(gpmodel.induce)==0
%     [M, S, V] = gp0_kuka_planar_dyn(gpmodel, m, s); return; end

D = length(m); E = size(gpmodel.targets,2);
%% Initialization for robot dynamics
persistent dynamics OPTIONS ctrlfcn u0 par jointlist njoint;
if isempty(dynamics)
    dynamics    = @dynamics_kp_nop;
    OPTIONS     = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
    ctrlfcn     = str2func('zoh');   
    par.dt = gpmodel.stepsize; par.delay = 0; par.tau = gpmodel.stepsize;
    jointlist   = gpmodel.jointi;
    njoint      = length(jointlist);
    u0          = cell(1,njoint);
end

%% Robot Dynamics
% n_span  = gpmodel.n_span;
q       = m(jointlist);
qdot    = m(jointlist + njoint);
tau     = m(end-njoint+1:end);

% for j = 1:njoint, u0{j} = @(t)ctrlfcn(tau(j,:),t,par); end
% [~, y] = ode45(dynamics, linspace(0,gpmodel.stepsize,n_span+1), m(1:(2*njoint)), OPTIONS, u0{:});

qddot   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,q,qdot,tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
[dqddotdq, dqddotdqdot, dqddotdtau, dqddotdqdq, dqddotdqdqdot, dqddotdqdtau, dqddotdqdotdqdot, dqddotdqdotdtau, dqddotdtaudtau] = ...
solveForwardDynamicsSecondDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,q,qdot,qddot,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
% for i = 2:n_span
%     q_tmp       = y(i,jointlist)';
%     qdot_tmp    = y(i,jointlist + njoint)';
%     qddot_tmp   = solveForwardDynamics(gpmodel.robot.A,gpmodel.robot.M,q_tmp,qdot_tmp,tau,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
%     [dqddotdq_temp, dqddotdqdot_temp, dqddotdtau_temp, dqddotdqdq_temp, dqddotdqdqdot_temp, dqddotdqdtau_temp, dqddotdqdotdqdot_temp, dqddotdqdotdtau_temp, dqddotdtaudtau_temp] = ...
%         solveForwardDynamicsSecondDerivatives_pilco(gpmodel.robot.A,gpmodel.robot.M,q_tmp,qdot_tmp,qddot_tmp,gpmodel.robot.G,gpmodel.Vdot0, gpmodel.robot.F);
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
% dqddotdtau          = dqddotdtau * gpmodel.stepsize/n_span;
dqddotdq            = dqddotdq * gpmodel.stepsize;
dqddotdqdot         = dqddotdqdot * gpmodel.stepsize;
dqddotdtau          = dqddotdtau * gpmodel.stepsize;
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
M(jointlist)            = qdot * gpmodel.stepsize;
M(jointlist + njoint)   = qddot * gpmodel.stepsize;

A                                        = zeros(E,D);
A(jointlist,jointlist + njoint)          = eye(njoint,njoint) * gpmodel.stepsize;
A(jointlist + njoint,jointlist)          = dqddotdq;
A(jointlist + njoint,jointlist + njoint) = dqddotdqdot;
A(jointlist + njoint,end-njoint+1:end)   = dqddotdtau;

dAdm = zeros(E,D,D);
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i)                = dFDdqdq(:,:,i);
    dAdm(jointlist + njoint,jointlist + njoint,i)       = dFDdqdotdq(:,:,i);
    dAdm(jointlist + njoint,end-njoint+1:end,i)         = dFDdtaudq(:,:,i);   
end
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i+njoint)             = dFDdqdqdot(:,:,i);
    dAdm(jointlist + njoint,jointlist + njoint,i+njoint)    = dFDdqdotdqdot(:,:,i);
    dAdm(jointlist + njoint,end-njoint+1:end,i+njoint)      = dFDdtaudqdot(:,:,i);   
end
for i = 1:njoint
    dAdm(jointlist + njoint,jointlist,i+end-njoint)             = dFDdqdtau(:,:,i);
    dAdm(jointlist + njoint,jointlist + njoint,i+end-njoint)    = dFDdqdotdtau(:,:,i);
    dAdm(jointlist + njoint,end-njoint+1:end,i+end-njoint)      = dFDdtaudtau(:,:,i);   
end


S =  A * s * (A');
S = (S + S')/2;

V = s * A';

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
            dSDds(:,:,i,j)    = A(:,i) * (A(:,j)') + A(:,j) * (A(:,i)');
            dSDds(:,:,j,i)    = dSDds(:,:,i,j);
            tempmatdVD(i,:)   = A(:,j)';
            tempmatdVD(j,:)   = tempmatdVD(j,:) + A(:,i)';
            dVDds(:,:,i,j)    = tempmatdVD;
            dVDds(:,:,j,i)    = tempmatdVD;
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
% dSdm = (dSDdm + permute(dSDdm, [2,1,3]))/2;
% dSds = (dSDds + permute(dSDds,[2,1,3,4]))/2;
% dSds = (dSds + permute(dSds,[1,2,4,3]))/2;
dVdm = permute(dAdm,[2,1,3]);
for i = 1:D
   dVdm(:,:,i) = s * dVdm(:,:,i); 
end
% dVds = (dVDds + permute(dVDds,[1,2,4,3]))/2;
%%

% % 5) vectorize derivatives
% dMds = reshape(dMds,[E D*D]);
% dSds = kron(A,A); % kron(A_g,A_g) == dSDds
% dSdm = reshape(dSdm,[E*E D]);
% dVds = reshape(dVds,[D*E D*D]);
% dVdm = reshape(dVdm,[D*E D]);
% 
% % 6) symmetrize variables dependent to variance
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