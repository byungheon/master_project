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
qdot0 = zeros(7,1);
qddot0 = zeros(7,1);

dt = 0.05;
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
trajectory = zeros(7,n_step+1);
trajectory(:,1) = x0s(1:7);

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
%     U = U(1:6) + kp * (q0(1:6) - x0s(1:6)) + kd * (-x0s(8:13));
%     U = [U;0];
    U(1:6) = U(1:6) + kp * (q0(1:6) - q(1:6)) + kd * (-dq(1:6));
    disp(U')
    U(7) = 0;
    ddq = solveForwardDynamics(A,M,q,dq,U,G, Vdot_0, robot.F);
    
    dq = dq + ddq * dt;
    q  = q + dq * dt;
    
    trajectory(:,i+1) = q;
end

 appVisualizeKUKA(trajectory);

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