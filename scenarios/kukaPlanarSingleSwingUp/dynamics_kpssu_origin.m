%% dynamics_kpssu_origin.m
% *Summary:* Implements ths ODE for simulating the cart-double pendulum 
% dynamics. 
%
%    function dz = dynamics_kpssu(t,z,f)
%
%
% *Input arguments:*
%
%		t     current time step (called from ODE solver)
%   z     state                                                    [6 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%   
%   dz    if 3 input arguments:      state derivative wrt time
%
%   Note: It is assumed that the state variables are of the following order:
%         q(3):        [rad]     angle of joints
%         dq(3):          [rad/s]   angular velocity of joints

function dz = dynamics_kpssu_origin(t,z,f1,f2)
%% Code

% set up the system
persistent robot_simul;
if(~isfield(robot_simul, 'A'))
   disp("robot initial construction");
   robot_simul = makeKukaR820_planar();
end
Vdot_0 = zeros(6,1); Vdot_0(6) = 9.82;
tau = zeros(3,1);
tau(1) = f1(t); 
tau(2) = f2(t);


q = z(1:3);
dq = z(4:6);
dz = zeros(6,1);
dz(1:3) = dq;
dz(4:6) = solveForwardDynamics(robot_simul.A,robot_simul.M,q,dq,tau,robot_simul.G, Vdot_0, robot_simul.F);
% disp(['odez:' num2str(dz(8:14)')]);
end