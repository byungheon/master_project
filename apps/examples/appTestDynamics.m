close all
clear
clc

%%
robot = makeKukaR820;

n = robot.dof;
m = 3;
A = robot.A;
M = robot.M;
G = robot.G;
Phi = robot.Phi;

q = rand(7,1);
qdot = rand(7,1);
qddot = rand(7,1);
tau   = rand(7,1);
%% Inverse Dynamics
[tau_from_inverse_dynamics, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G);


%% Forward Dynamics
disp(' ')

% tic
% for i = 1:1000
% [tau_from_inverse_dynamics, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G);
% end
% toc
vdot_0 = [0;0;0;0;0;9.81];
q = rand(7,1);
qdot = rand(7,1);
qddot = rand(7,1);
tic
for i =1:1000
qddot_f_1 = solveForwardDynamics(A,M,q,qdot,tau(:,1),G, vdot_0);
end
toc

tic
for i = 1:1000
qddot_f_2 = solveForwardDynamics2(A,M,q,qdot,tau(:,1),G, vdot_0);
end
toc
temp = [qddot_f_1 qddot_f_2];
qddot_f_1 - qddot_f_2
%%
close all;
angle = [0;pi/8;0;-3*pi/8;pi/2;pi/2;0];
appVisualizeKUKA(angle);