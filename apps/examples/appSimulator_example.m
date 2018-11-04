close all
clear
clc

%%
robot = makeKukaR820;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
Phi = robot.Phi;

q0 = [0;pi/8;0;-3*pi/8;pi/2;pi/2;pi/2];
qdot0 = zeros(n,1);
qddot0 = zeros(n,1);

dt = 0.01;
T  = 10;
n_step = ceil(T/dt);
dynamics = @dynamics_kssu;
OPTIONS = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
ctrlfcn = str2func('zoh');
nU = 6;
u0 = cell(1,nU);
x0s = [q0;qdot0];
par.dt = dt; par.delay = 0; par.tau = dt;
Vdot_0 = zeros(6,1);
Vdot_0(6) = 9.82;
trajectory1 = zeros(7,n_step+1);
trajectory2 = zeros(7,n_step+1);

trajectory1(:,1) = x0s(1:7);
trajectory2(:,1) = x0s(1:7);

% for i = 1:n_step
%     U = solveInverseDynamics(A, M, x0s(1:7), zeros(7,1), zeros(7,1), G, Vdot_0);
%     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
%     for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
%     [T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     x0s = y(3,:)';
%     trajectory(:,i+1) = x0s(1:7);
% end

kp = 2;
kd = 0.1;
q  = x0s(1:7);
dq = x0s(8:14);
for i = 1:n_step
    U = solveInverseDynamics(A, M, q, zeros(7,1), zeros(7,1), G, Vdot_0);
    disp(U');
% %     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
% %     U = [U;0];
%     U(1:6) = U(1:6) + kp * (q0(1:6) - q(1:6)) + kd * (-dq(1:6));
%     disp(U')
    U(7) = 0;
    ddq = solveForwardDynamics(A,M,q,dq,zeros(7,1),G, Vdot_0, robot.F);     
    dq = dq + ddq * dt;
    q  = q + dq * dt;
    trajectory1(:,i+1) = q;
    
    U = zeros(6,1);
    for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
    [T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     ddq = solveForwardDynamics(A,M,x0s(1:7),x0s(8:14),zeros(7,1),G, Vdot_0, robot.F);
    disp('---------------');
%     dqstep      = (y(3,1:7)' - x0s(1:7))./x0s(8:14); 
%     dqdotstep   = (y(3,8:14)' - x0s(8:14))./ddq;
%     disp(dqstep');
%     disp(dqdotstep');
% %     disp(ddq');
% %     disp(-(x0s(8:14)'-y(3,8:14))/dt);
%     disp('------------------------------------------------------------------------------------');
    x0s = y(3,:)';
    trajectory2(:,i+1) = x0s(1:7);
end
% 
 appVisualizeKUKA(trajectory1);
%  appVisualizeKUKA(trajectory2);

%%
%%
clear all
close all
clc
robot = makeKukaR820_prior;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
Phi = robot.Phi;

q0 = [0;pi/8;0;-3*pi/8;pi/2;pi/2];
qdot0 = zeros(n,1);
qddot0 = zeros(n,1);

dt = 0.001;
T  = 5;
n_step = ceil(T/dt);
dynamics = @dynamics_kuka_6dof;
OPTIONS = odeset('RelTol', 1e-5, 'AbsTol', 1e-5);
ctrlfcn = str2func('zoh');
nU = 6;
u0 = cell(1,nU);
x0s = [q0;qdot0];
par.dt = dt; par.delay = 0; par.tau = dt;
Vdot_0 = zeros(6,1);
Vdot_0(6) = 9.82;
trajectory1 = zeros(6,n_step+1);
trajectory2 = zeros(6,n_step+1);

trajectory1(:,1) = x0s(1:6);
trajectory2(:,1) = x0s(1:6);


kp = 2;
kd = 0.1;
q  = x0s(1:6);
dq = x0s(7:12);
for i = 1:n_step
%     U = solveInverseDynamics(A, M, q, zeros(7,1), zeros(7,1), G, Vdot_0);
% %     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
% %     U = [U;0];
%     U(1:6) = U(1:6) + kp * (q0(1:6) - q(1:6)) + kd * (-dq(1:6));
%     disp(U')
    
    ddq = solveForwardDynamics(A,M,q,dq,zeros(6,1),G, Vdot_0, robot.F);     
    dq = dq + ddq * dt;
    q  = q + dq * dt;
    trajectory1(:,i+1) = q;
%     disp(dq');
    
    U = zeros(6,1);
    for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
    [T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     ddq = solveForwardDynamics(A,M,x0s(1:7),x0s(8:14),zeros(7,1),G, Vdot_0, robot.F);
%     disp('---------------');
%     dqstep      = (y(3,1:7)' - x0s(1:7))./x0s(8:14); 
%     dqdotstep   = (y(3,8:14)' - x0s(8:14))./ddq;
%     disp(dqstep');
%     disp(dqdotstep');
% %     disp(ddq');
% %     disp(-(x0s(8:14)'-y(3,8:14))/dt);
%     disp('------------------------------------------------------------------------------------');
    x0s = y(3,:)';
    trajectory2(:,i+1) = x0s(1:6);
end
% 
%  appVisualizeKUKA_6dof(trajectory1);
%  appVisualizeKUKA(trajectory2);
for i = 1:6
figure(i);
plot(trajectory1(i,:));
hold on;
plot(trajectory2(i,:));
legend('simple','ode');
end

%%
clear all
close all
clc
robot = makeKukaR820_planar;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;

q0 = zeros(n,1);
qdot0 = zeros(n,1);
qddot0 = zeros(n,1);

dt = 0.01;
T  = 5;
n_step = ceil(T/dt);
dynamics1 = @dynamics_kpssu;
dynamics2 = @dynamics_kpssu_origin;
OPTIONS = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
ctrlfcn = str2func('zoh');
nU = 2;
u0 = cell(1,nU);
x0s1 = [q0;qdot0];
x0s2 = [q0;qdot0];
par.dt = dt; par.delay = 0; par.tau = dt;
Vdot_0 = zeros(6,1);
Vdot_0(6) = 9.82;
trajectory1 = zeros(3,n_step+1);
trajectory2 = zeros(3,n_step+1);

trajectory1(:,1) = x0s1(1:3);
trajectory2(:,1) = x0s2(1:3);

% for i = 1:n_step
%     U = solveInverseDynamics(A, M, x0s(1:7), zeros(7,1), zeros(7,1), G, Vdot_0);
%     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
%     for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
%     [T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     x0s = y(3,:)';
%     trajectory(:,i+1) = x0s(1:7);
% end

kp = 2;
kd = 0.1;
q  = x0s1(1:3);
dq = x0s1(4:6);
for i = 1:n_step
%     U = solveInverseDynamics(A, M, q, zeros(3,1), zeros(3,1), G, Vdot_0);
%     disp(U');
% %     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
% %     U = [U;0];
%     U(1:6) = U(1:6) + kp * (q0(1:6) - q(1:6)) + kd * (-dq(1:6));
%     disp(U')
    
%     ddq = solveForwardDynamics(A,M,q,dq,zeros(3,1),G, Vdot_0, robot.F);     
%     dq = dq + ddq * dt;
%     q  = q + dq * dt;
    
    U = zeros(2,1);
    for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
    [T1 y1] = ode45(dynamics1, [0 dt/2 dt], x0s1, OPTIONS, u0{:});
    x0s1 = y1(3,:)';
    [T2 y2] = ode45(dynamics2, [0 dt/2 dt], x0s2, OPTIONS, u0{:});
    x0s2 = y2(3,:)';

    tempq = x0s1(1:3);
    tempq(3) = tempq(3) - tempq(1) + tempq(2);
    trajectory1(:,i+1) = tempq;
    trajectory2(:,i+1) = x0s2(1:3);
    
end
close all
for i = 1:3
figure(i);
plot(trajectory1(i,:)-0.2);
hold on;
plot(trajectory2(i,:));
legend('simple','ode');
end
% 
%  appVisualizeKUKA_planar(trajectory1);
%  appVisualizeKUKA_planar(trajectory2);

%%
clear all
close all
clc
robot = makeKukaR820_planar_prior;
rand('seed',10);
n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;

q0 = zeros(n,1);
qdot0 = zeros(n,1);
qddot0 = zeros(n,1);

dt = 0.005;
T  = 5;
n_step = ceil(T/dt);
dynamics1 = @dynamics_kp_nop;
OPTIONS = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
ctrlfcn = str2func('zoh');
nU = 2;
u0 = cell(1,nU);
x0s1 = [q0;qdot0];
x0s2 = [q0;qdot0];
par.dt = dt; par.delay = 0; par.tau = dt;
Vdot_0 = zeros(6,1);
Vdot_0(6) = 9.82;
trajectory1 = zeros(2*n,n_step+1);
trajectory2 = zeros(2*n,n_step+1);

trajectory1(:,1) = x0s1;
trajectory2(:,1) = x0s2;

kp = 2;
kd = 0.1;
q  = x0s1(1:n);
dq = x0s1(n+1:2*n);
x0s2 = x0s1;
tic
for i = 1:n_step
%     U = solveInverseDynamics(A, M, q, zeros(3,1), zeros(3,1), G, Vdot_0);
%     disp(U');
% %     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
% %     U = [U;0];
%     U(1:6) = U(1:6) + kp * (q0(1:6) - q(1:6)) + kd * (-dq(1:6));
%     disp(U')
%     U = zeros(2,1);
%     U = -[60;70] + [120;140].* rand(2,1);
%     ddq = solveForwardDynamics(A,M,q,dq,U,G, Vdot_0, robot.F);     
%     dq = dq + ddq * dt;
%     q  = q + dq * dt;
%     dq = dq + ddq * dt;
%     dq = dq + ddq * dt;
    
    U = -[60;70] + [120;140].* rand(2,1);
    
    for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
    [T1 y1] = ode45(dynamics1, [0 dt/2 dt], x0s1, OPTIONS, u0{:});
    x0s1 = y1(3,:)';
    
    [T2 y2] = ode45(@(t,y)dynamics_kp_nop_not(t,y,U(1),U(2)), [0 dt/2 dt], x0s2, OPTIONS);
    x0s2 = y2(3,:)';
    
%     if y1(3,:) == y2(3,:)
%        disp('yeah!!');
%     else
%         disp('fuck');
%     end
%     trajectory1(:,i+1) = x0s1;
%     trajectory2(:,i+1) = [q;dq];
    
end
toc
tic
for i = 1:n_step

    U = -[60;70] + [120;140].* rand(2,1);

    [T2 y2] = ode45(@(t,y)dynamics_kp_nop_not(t,y,U(1),U(2)), [0 dt/2 dt], x0s2, OPTIONS);
    x0s2 = y2(3,:)';
    
%     
%     trajectory1(:,i+1) = x0s1;
%     trajectory2(:,i+1) = [q;dq];
    
end
toc
% close all
% for i = 1:2*n
% figure(i);
% plot(trajectory1(i,:));
% hold on;
% plot(trajectory2(i,:));
% legend('ode','simple');
% end
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