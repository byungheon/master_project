%% dynamics_kssu.m
% *Summary:* Implements ths ODE for simulating the cart-double pendulum 
% dynamics. 
%
%    function dz = dynamics_kssu(t,z,f)
%
%
% *Input arguments:*
%
%		t     current time step (called from ODE solver)
%   z     state                                                    [14 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%   
%   dz    if 3 input arguments:      state derivative wrt time
%
%   Note: It is assumed that the state variables are of the following order:
%         q(7):        [rad]     angle of joints
%         dq(7):          [rad/s]   angular velocity of joints

function dz = dynamics_kssu(t,z,f1,f2,f3,f4,f5,f6)
%% Code

% set up the system
persistent robot_simul;
if(~isfield(robot_simul, 'A'))
   disp("robot initial construction");
   robot_simul = makeKukaR820;
end
Vdot_0 = zeros(6,1); Vdot_0(6) = 9.82;
tau = zeros(7,1);
tau(1) = f1(t); 
tau(2) = f2(t);
tau(3) = f3(t);
tau(4) = f4(t);
tau(5) = f5(t);
tau(6) = f6(t);

q = z(1:7);
dq = z(8:14);
dz = zeros(14,1);
dz(1:7) = dq;
dz(8:14) = solveForwardDynamics(robot_simul.A,robot_simul.M,q,dq,tau(:),robot_simul.G, Vdot_0, robot_simul.F, robot_simul.Sigmoid);
% disp(['odez:' num2str(dz(8:14)')]);
end