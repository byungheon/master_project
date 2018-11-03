clear all;
close all
clc;
settings_NOGP;
dynmodel.targets = zeros(3,6);

robot = makeKukaR820_planar;

% n = robot.dof;
% A = robot.A;
% M = robot.M;
% G = robot.G;
% Phi = robot.Phi;

M = zeros(14,1);
S = zeros(14,14);
m0 = [0;0;0;0;0;0];
q0 = m0;
% qdot0 = zeros(n,1);
% qddot0 = zeros(n,1);
% dt = dynmodel.dt;
T  = 5;
n_step = ceil(T/dt);

s0 = diag([0.01;0.01;0.01;0.01;0.01;0.01].^2);
tau0 = [0;0];
taus0 = diag([0.01;0.01].^2);
cell_M = zeros(6,n_step+1);
cell_S = zeros(6,6,n_step+1);
cell_M(:,1) = m0;
cell_S(:,:,1) = s0;

dynamics = @dynamics_NOGP_learning;
OPTIONS = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
ctrlfcn = str2func('zoh');
nU = 2;
u0 = cell(1,nU);
x0s = q0;
par.dt = dt; par.delay = 0; par.tau = dt;

trajectory1 = zeros(6,n_step+1);
trajectory1(:,1) = x0s;
for step_i = 1:n_step
M(1:6) = m0;
S(1:6,1:6) = s0;
M(7:8) = tau0;
S(7:8,7:8) = taus0;
i = 1:8;
k = 9:14;
j = 1:8;
[M(k), S(k,k), C] = gp1_NOGP(dynmodel, M(i), S(i,i));
q = S(j,i)*C; S(j,k) = q; S(k,j) = q';

% 4) Compute distribution of the next state -----------------------------------
P = [zeros(6,8) eye(6)]; P(1:6,1:6) = eye(6);
Mnext = P*M; Snext = P*S*P'; Snext = (Snext+Snext')/2;

m0 = Mnext;
s0 = Snext;
cell_M(:,step_i+1)    = Mnext;
cell_S(:,:,step_i+1)  = Snext;

U = zeros(2,1);
for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
[T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     ddq = solveForwardDynamics(A,M,x0s(1:7),x0s(8:14),zeros(7,1),G, Vdot_0, robot.F);

x0s = y(3,:)';
trajectory1(:,step_i+1) = x0s;

% disp(x0s');
% disp(Mnext');
% disp('------------------');

end






for i = 1:6
if ~ishandle(i); figure(i); else set(0,'CurrentFigure',i); end
  clf(i); 
  std = cell_S(i,i,:);
  std = std(:)';
  errorbar(0:n_step,cell_M(i,:),2*sqrt(std)); drawnow;
%   plot(0:n_step,cell_M(i,:),'b');
    hold on;
  plot(0:n_step,trajectory1(i,:),'r');
end
%%
clear all;
close all
clc;

settings_kpssu;

robot = makeKukaR820_planar_prior;
n = robot.dof;
dynmodel.targets = zeros(3,2*n);

% n = robot.dof;
% A = robot.A;
% M = robot.M;
% G = robot.G;
% Phi = robot.Phi;

M = zeros(4*n+2,1);
S = zeros(4*n+2,4*n+2);
m0 = zeros(2*n,1);
q0 = m0;
% qdot0 = zeros(n,1);
% qddot0 = zeros(n,1);
% dt = dynmodel.dt;
T  = 5;
n_step = ceil(T/dt);

s0 = diag(0.01^2 * ones(2*n,1));
tau0 = [0;0];
taus0 = diag([0.01;0.01].^2);
cell_M = zeros(2*n,n_step+1);
cell_S = zeros(2*n,2*n,n_step+1);
cell_M(:,1) = m0;
cell_S(:,:,1) = s0;

dynamics = @dynamics_kp_nop;
OPTIONS = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
ctrlfcn = str2func('zoh');
nU = 2;
u0 = cell(1,nU);
x0s = q0;
par.dt = dt; par.delay = 0; par.tau = dt;

trajectory1 = zeros(2*n,n_step+1);
trajectory1(:,1) = x0s;
for step_i = 1:n_step
M(1:2*n) = m0;
S(1:2*n,1:2*n) = s0;
M(2*n+1:2*n+2) = tau0;
S(2*n+1:2*n+2,2*n+1:2*n+2) = taus0;
i = 1:2*n+2;
k = 2*n+2+1:4*n+2;
j = 1:2*n+2;
[M(k), S(k,k), C] = gp1_test(dynmodel, M(i), S(i,i));
q = S(j,i)*C; S(j,k) = q; S(k,j) = q';

% 4) Compute distribution of the next state -----------------------------------
P = [zeros(2*n,2*n+2) eye(2*n)]; P(1:2*n,1:2*n) = eye(2*n);
Mnext = P*M; Snext = P*S*P'; Snext = (Snext+Snext')/2;

m0 = Mnext;
s0 = Snext;
cell_M(:,step_i+1)    = Mnext;
cell_S(:,:,step_i+1)  = Snext;

U = zeros(2,1);
for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
[T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     ddq = solveForwardDynamics(A,M,x0s(1:7),x0s(8:14),zeros(7,1),G, Vdot_0, robot.F);

x0s = y(3,:)';
trajectory1(:,step_i+1) = x0s;

% disp(x0s');
% disp(Mnext');
% disp('------------------');

end

for i = 1:2*n
if ~ishandle(i); figure(i); else set(0,'CurrentFigure',i); end
  clf(i); 
  std = cell_S(i,i,:);
  std = std(:)';
  errorbar(0:n_step,cell_M(i,:),2*sqrt(std)); drawnow;
%   plot(0:n_step,cell_M(i,:),'b');
    hold on;
  plot(0:n_step,trajectory1(i,:),'r');
end

%%
%%
clear all;
close all
clc;

settings_kpssu;
% load('C:\Users\my\Desktop\bh\master_project\data\KukaPlanarSingleSwingUp_Temp_PILCO_2_H250.mat');

robot = makeKukaR820_planar_prior;
n = robot.dof;

% n = robot.dof;
% A = robot.A;
% M = robot.M;
% G = robot.G;
% Phi = robot.Phi;

m0 = - 1 + 2 * rand(2*n,1);
s0 = diag(0.1^2 * ones(2*n,1));
tau0 = - 10 + 20 * rand(2,1);
taus0 = diag([0.1;0.1].^2);

M0 = [m0;tau0];
S0 = diag(0.01^2 * ones(2*n+2,1));

dt = 0.00001;
% dynmodel.stepsize = dt;
% x_start = randn();
x_start = 0.0665;
x       = x_start:dt:(x_start+10*dt);

D = 2*n+2;
E = 2*n;
% D = 11;
% E = 6;
M0 = -1 + 2 * rand(D,1);
% M0 = [    0.4540
%     0.3113
%    -0.8053
%     0.7822
%     0.0252
%    -0.6122
%    -0.9575
%     0.7448
%    -0.2473
%    -0.4991
%     0.1352];
S0 = diag(0.01^2 * ones(D,1));

M = zeros(E,length(x));
S = zeros(E,E,length(x));
V = zeros(D,E,length(x));
dMdm = zeros(E,D,length(x));
dSdm = zeros(E,E,D,length(x));
dVdm = zeros(D,E,D,length(x));
dMds = zeros(E,D,D,length(x));
dSds = zeros(E,E,D,D,length(x));
dVds = zeros(D,E,D,D,length(x));
m_i = 2;

% dynmodel.inputs     = rand(100,D);
dynmodel.targets    = rand(100,E);
% dynmodel.hyp        = rand(D+2,E);
% 
dynmodel.n_span = 200;
for m_i = 1:D
    for i = 1:length(x)
    m = M0;
    s = S0;
    
    m(m_i) = m(m_i) + x(i);

    [M(:,i), S(:,:,i), V(:,:,i), dMdm(:,:,i), dSdm(:,:,:,i), dVdm(:,:,:,i), dMds(:,:,:,i), dSds(:,:,:,:,i), dVds(:,:,:,:,i)] = gp1d_test(dynmodel, m, s);

    end

%dMdm    
    numerical = zeros(E,length(x));
    for i = 1:E
       gradient_tmp = M(i,:);
       gradient_tmp = gradient_tmp(:)';
       numerical(i,:) = gradient(gradient_tmp,dt);
    end
    numerical = numerical(:,2:end-1);
    analytic  = dMdm(:, m_i, :);
    analytic  = analytic(:,2:end-1);
    numerical
    disp('--------------------------');
    analytic
    disp('--------------------------');
    analytic - numerical 
    disp('--------------------------');
    
%dSdm
%     numerical = zeros(E*E,length(x));
%     for i = 1:E
%         for j=1:E
%            gradient_tmp = S(i,j,:);
%            gradient_tmp = gradient_tmp(:)';
%            numerical((j-1)*E+i,:) = gradient(gradient_tmp,dt);
%         end
%     end
%     numerical = numerical(:,2:end-1);
%     analytic  = dSdm(:, :, m_i, :);
%     analytic  = reshape(analytic,[E*E,size(analytic,4)]);
%     analytic  = analytic(:,2:end-1);
% %     numerical
% %     disp('--------------------------');
%     analytic - numerical 
%     disp('--------------------------');

%dVdm
%     numerical = zeros(D*E,length(x));
%     for i = 1:D
%         for j=1:E
%            gradient_tmp = V(i,j,:);
%            gradient_tmp = gradient_tmp(:)';
%            numerical((j-1)*D+i,:) = gradient(gradient_tmp,dt);
%         end
%     end
%     numerical = numerical(:,2:end-1);
%     analytic  = dVdm(:, :, m_i, :);
%     analytic  = reshape(analytic,[D*E,size(analytic,4)]);
%     analytic  = analytic(:,2:end-1);numerical
%     disp('--------------------------');
%     analytic - numerical 
%     disp('--------------------------');
    disp('--------------------------');
end

% tic
% error = zeros(E,length(x)-2);
% for sc_i = 1:D
%     for sr_i = 1:D
%         for i = 1:length(x)
%         m = M0;
%         s = S0;
% 
%     %     m(m_i) = m(m_i) + x(i);
%         s(sc_i,sr_i) = s(sc_i,sr_i) + x(i);
%         s(sr_i,sc_i) = s(sc_i,sr_i);
%         
%         [M(:,i), S(:,:,i), V(:,:,i), dMdm(:,:,i), dSdm(:,:,:,i), dVdm(:,:,:,i), dMds(:,:,:,i), dSds(:,:,:,:,i), dVds(:,:,:,:,i)] = gp1d_odetest(dynmodel, m, s);
%         
%         end
% %         dMds
% %         numerical = zeros(E,length(x));
% %         for i = 1:E
% %            gradient_tmp = M(i,:);
% %            gradient_tmp = gradient_tmp(:)';
% %            numerical(i,:) = gradient(gradient_tmp,dt);
% %         end
% %         numerical = numerical(:,2:end-1);
% %         if sc_i ~= sr_i
% %             analytic  = dMds(:, sc_i, sr_i, :);
% %         else
% %             analytic  = dMds(:, sc_i, sr_i, :);
% %         end
% %         
% %         analytic  = reshape(analytic,[E,size(analytic,4)]);
% %         analytic  = analytic(:,2:end-1);numerical
% %     disp('--------------------------');
% %     analytic - numerical 
% %     disp('--------------------------');
% %         if ~isreal(analytic - numerical)
% %             if ~isreal(analytic)
% %                 disp('analytic complex!!');
% %             else
% %                 disp('numerical complex!!');
% %             end
% %             keyboard;
% %         end
% %         error = error + abs(analytic - numerical);
%     
% %dSds
% %         numerical = zeros(E*E,length(x));
% %         for i = 1:E
% %             for j=1:E
% %                gradient_tmp = S(i,j,:);
% %                gradient_tmp = gradient_tmp(:)';
% %                numerical((j-1)*E+i,:) = gradient(gradient_tmp,dt);
% %             end
% %         end
% %         numerical = numerical(:,2:end-1);
% %         analytic  = dSds(:, :, sc_i, sr_i, :);
% %         analytic  = reshape(analytic,[E*E,size(analytic,5)]);
% %         analytic  = analytic(:,2:end-1);
% % %         numerical
% % %     disp('--------------------------');
% %     analytic - numerical 
% %     disp('--------------------------');
% 
% % %dVds
% %         numerical = zeros(D*E, length(x));
% %         for i = 1:D
% %             for j = 1:E
% %                gradient_tmp = V(i,j,:);
% %                gradient_tmp = gradient_tmp(:)';
% %                numerical((j-1)*D+i,:) = gradient(gradient_tmp,dt);
% %             end
% %         end
% %         numerical = numerical(:,2:end-1);
% %         analytic  = dVds(:, :, sc_i,sr_i, :);
% %        
% %         
% %         analytic  = reshape(analytic,[D*E,size(analytic,5)]);
% %         analytic  = analytic(:,2:end-1);
% % %         numerical
% % %     disp('--------------------------');
% %     analytic - numerical 
% %     disp('--------------------------');
% %         disp('--------------------------');
%     end
% end
% toc
%%
%%
clear all;
close all
clc;

settings_NOGP;
% load('C:\Users\my\Desktop\bh\master_project\data\KukaPlanarSingleSwingUp_Temp_PILCO_2_H250.mat');

robot = makeKukaR820_planar;
n = robot.dof;
n_control = 2;
% n = robot.dof;
% A = robot.A;
% M = robot.M;
% G = robot.G;
% Phi = robot.Phi;

m0 = - 1 + 2 * rand(2*n,1);
s0 = diag(0.1^2 * ones(2*n,1));
tau0 = - 10 + 20 * rand(2,1);
taus0 = diag(0.1^2 * ones(n_control,1));


dt = 0.00001;
% dynmodel.stepsize = dt;
% x_start = randn();
x_start = 0.0665;
x       = x_start:dt:(x_start+10*dt);

D = 2*n+n_control;
E = 2*n;
% D = 11;
% E = 6;
M0 = -1 + 2 * rand(D,1);
% M0 = [    0.4540
%     0.3113
%    -0.8053
%     0.7822
%     0.0252
%    -0.6122
%    -0.9575
%     0.7448
%    -0.2473
%    -0.4991
%     0.1352];
S0 = diag(0.01^2 * ones(D,1));

M = zeros(E,length(x));
S = zeros(E,E,length(x));
V = zeros(D,E,length(x));
dMdm = zeros(E,D,length(x));
dSdm = zeros(E,E,D,length(x));
dVdm = zeros(D,E,D,length(x));
dMds = zeros(E,D,D,length(x));
dSds = zeros(E,E,D,D,length(x));
dVds = zeros(D,E,D,D,length(x));

aM = zeros(E,length(x));
aS = zeros(E,E,length(x));
aV = zeros(D,E,length(x));
adMdm = zeros(E,D,length(x));
adSdm = zeros(E,E,D,length(x));
adVdm = zeros(D,E,D,length(x));
adMds = zeros(E,D,D,length(x));
adSds = zeros(E,E,D,D,length(x));
adVds = zeros(D,E,D,D,length(x));

% dynmodel.inputs     = rand(100,D);
dynmodel.targets    = rand(100,E);
% dynmodel.hyp        = rand(D+2,E);
% 
dynmodel.n_span = 200;
% for m_i = 1:D
%     for i = 1:length(x)
%     m = M0;
%     s = S0;
%     
%     m(m_i) = m(m_i) + x(i);
%     
%     [M(:,i), S(:,:,i), V(:,:,i), dMdm(:,:,i), dSdm(:,:,:,i), dVdm(:,:,:,i), dMds(:,:,:,i), dSds(:,:,:,:,i), dVds(:,:,:,:,i)] = gp1d_NOGP_MINE(dynmodel, m, s);
% 
%     end
% 
% %dMdm    
% %     numerical = zeros(E,length(x));
% %     for i = 1:E
% %        gradient_tmp = M(i,:);
% %        gradient_tmp = gradient_tmp(:)';
% %        numerical(i,:) = gradient(gradient_tmp,dt);
% %     end
% %     numerical = numerical(:,2:end-1);
% %     analytic  = dMdm(:, m_i, :);
% %     analytic  = analytic(:,2:end-1);
% % %     numerical
% % %     disp('--------------------------');
% % %     analytic
% % %     disp('--------------------------');
% %     analytic - numerical 
% %     disp('--------------------------');
%     
% %dSdm
% %     numerical = zeros(E*E,length(x));
% %     for i = 1:E
% %         for j=1:E
% %            gradient_tmp = S(i,j,:);
% %            gradient_tmp = gradient_tmp(:)';
% %            numerical((j-1)*E+i,:) = gradient(gradient_tmp,dt);
% %         end
% %     end
% %     numerical = numerical(:,2:end-1);
% %     analytic  = dSdm(:, :, m_i, :);
% %     analytic  = reshape(analytic,[E*E,size(analytic,4)]);
% %     analytic  = analytic(:,2:end-1);
% %     numerical
% %     disp('--------------------------');
% %     analytic - numerical 
% %     disp('--------------------------');
% 
% %dVdm
%     numerical = zeros(D*E,length(x));
%     for i = 1:D
%         for j=1:E
%            gradient_tmp = V(i,j,:);
%            gradient_tmp = gradient_tmp(:)';
%            numerical((j-1)*D+i,:) = gradient(gradient_tmp,dt);
%         end
%     end
%     numerical = numerical(:,2:end-1);
%     analytic  = dVdm(:, :, m_i, :);
%     analytic  = reshape(analytic,[D*E,size(analytic,4)]);
%     analytic  = analytic(:,2:end-1);numerical
%     disp('--------------------------');
%     analytic - numerical 
%     disp('--------------------------');
%     disp('--------------------------');
% end

tic
error = zeros(E,length(x)-2);
for sc_i = 1:D
    for sr_i = 1:D
        for i = 1:length(x)
        m = M0;
        s = S0;

        s(sc_i,sr_i) = s(sc_i,sr_i) + x(i);
%         s(sr_i,sc_i) = s(sc_i,sr_i);
        
        [M(:,i), S(:,:,i), V(:,:,i), dMdm(:,:,i), dSdm(:,:,:,i), dVdm(:,:,:,i), dMds(:,:,:,i), dSds(:,:,:,:,i), dVds(:,:,:,:,i)] = gp1d_NOGP_MINE(dynmodel, m, s);
        
        end
        for i = 1:length(x)
        m = M0;
        s = S0;

%         s(sc_i,sr_i) = s(sc_i,sr_i) + x(i);
        s(sr_i,sc_i) = s(sr_i,sc_i) + x(i);
        
        [aM(:,i), aS(:,:,i), aV(:,:,i), adMdm(:,:,i), adSdm(:,:,:,i), adVdm(:,:,:,i), adMds(:,:,:,i), adSds(:,:,:,:,i), adVds(:,:,:,:,i)] = gp1d_NOGP_MINE(dynmodel, m, s);
        
        end
% %         dMds
%         numerical = zeros(E,length(x));
%         for i = 1:E
%            gradient_tmp = M(i,:);
%            gradient_tmp = gradient_tmp(:)';
%            numerical(i,:) = gradient(gradient_tmp,dt);
%         end
%         numerical = numerical(:,2:end-1);
%         if sc_i ~= sr_i
%             analytic  = dMds(:, sc_i, sr_i, :);
%         else
%             analytic  = dMds(:, sc_i, sr_i, :);
%         end
%         
%         analytic  = reshape(analytic,[E,size(analytic,4)]);
%         analytic  = analytic(:,2:end-1);
%         numerical
%     disp('--------------------------');
%     analytic - numerical 
%     disp('--------------------------');
%         if ~isreal(analytic - numerical)
%             if ~isreal(analytic)
%                 disp('analytic complex!!');
%             else
%                 disp('numerical complex!!');
%             end
%             keyboard;
%         end
%         error = error + abs(analytic - numerical);
    
%dSds
        numerical = zeros(E*E,length(x));
        for i = 1:E
            for j=1:E
               gradient_tmp = S(i,j,:);
               gradient_tmp = gradient_tmp(:)';
               numerical((j-1)*E+i,:) = gradient(gradient_tmp,dt);
            end
        end
        numerical = numerical(:,2:end-1);
        
        numerical2 = zeros(E*E,length(x));
        for i = 1:E
            for j=1:E
               gradient_tmp = aS(i,j,:);
               gradient_tmp = gradient_tmp(:)';
               numerical2((j-1)*E+i,:) = gradient(gradient_tmp,dt);
            end
        end
        numerical2 = numerical2(:,2:end-1);

        if sc_i ~= sr_i
            numerical = numerical + numerical2;
        end
        analytic  = dSds(:, :, sc_i, sr_i, :);
        analytic  = reshape(analytic,[E*E,size(analytic,5)]);
        analytic  = analytic(:,2:end-1);
%         numerical
%     disp('--------------------------');
    analytic - numerical 
    disp('--------------------------');

% %dVds
%         numerical = zeros(D*E, length(x));
%         for i = 1:D
%             for j = 1:E
%                gradient_tmp = V(i,j,:);
%                gradient_tmp = gradient_tmp(:)';
%                numerical((j-1)*D+i,:) = gradient(gradient_tmp,dt);
%             end
%         end
%         numerical = numerical(:,2:end-1);
%         analytic  = dVds(:, :, sc_i,sr_i, :);
%        
%         
%         analytic  = reshape(analytic,[D*E,size(analytic,5)]);
%         analytic  = analytic(:,2:end-1);
% %         numerical
% %     disp('--------------------------');
%     (analytic - numerical)'
% %     disp('--------------------------');
%         disp('--------------------------');
    end
end
toc

%%
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