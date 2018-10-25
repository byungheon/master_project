close all
clear
clc

%%
robot = makeTestRobot();

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
Phi = robot.Phi;

vdot_0 = [0;0;0;0;0;9.82];

q = [0.39 1.39 2.39 3.39 4.39 5.39]';
qdot = [0.38 0.38 0.38 0.38 0.38 0.38]';
qddot = [0.18 0.18 0.18 0.18 0.18 0.18]';

%% First Derivative Dynamics
tic
[tau, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G,vdot_0);


[dtaudq1, dtaudqdot2, dtaudqddot3] = solveInverseDynamicsDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,vdot_0)
toc
%% Second Derivative Dynamics
tic
[tau, T, V, Vdot, F] = solveInverseDynamics(A,M,q,qdot,qddot,G, vdot_0);


[dtaudq, dtaudqdot, dtaudqddot, dtaudqdq, dtaudqdqdot, dtaudqdqddot, dtaudqdotdqdot, dtaudqdotdqddot, dtaudqddotdqddot] = solveInverseDynamicsSecondDerivatives_pilco(A,M,q,qdot,G,T,V,Vdot,F,vdot_0);

toc
%%
tic
for i =1 :100
    [dqddotdq, dqddotdqdot, dqddotdtau, dqddotdqdq, dqddotdqdqdot, dqddotdqdtau, dqddotdqdotdqdot, dqddotdqdotdtau, dqddotdtaudtau] = solveForwardDynamicsSecondDerivatives_pilco(A,M,q,qdot,qddot,G,vdot_0);
end
toc
