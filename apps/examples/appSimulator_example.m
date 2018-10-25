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

dt = 0.01;
T  = 10;
n_step = ceil(T/dt);
dynamics = @dynamics_kuka_6dof;
OPTIONS = odeset('RelTol', 1e-3, 'AbsTol', 1e-3);
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
    disp(dq');
%     U = zeros(6,1);
%     for j = 1:nU, u0{j} = @(t)ctrlfcn(U(j,:),t,par); end
%     [T y] = ode45(dynamics, [0 dt/2 dt], x0s, OPTIONS, u0{:});
%     ddq = solveForwardDynamics(A,M,x0s(1:7),x0s(8:14),zeros(7,1),G, Vdot_0, robot.F);
%     disp('---------------');
%     dqstep      = (y(3,1:7)' - x0s(1:7))./x0s(8:14); 
%     dqdotstep   = (y(3,8:14)' - x0s(8:14))./ddq;
%     disp(dqstep');
%     disp(dqdotstep');
% %     disp(ddq');
% %     disp(-(x0s(8:14)'-y(3,8:14))/dt);
%     disp('------------------------------------------------------------------------------------');
%     x0s = y(3,:)';
%     trajectory2(:,i+1) = x0s(1:7);
end
% 
 appVisualizeKUKA_6dof(trajectory1);
%  appVisualizeKUKA(trajectory2);


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