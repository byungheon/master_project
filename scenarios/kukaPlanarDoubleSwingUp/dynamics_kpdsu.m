%% dynamics_kpssu.m
% *Summary:* Implements ths ODE for simulating the cart-double pendulum 
% dynamics. 
%
%    function dz = dynamics_kpssu(t,z,f)
%
%
% *Input arguments:*
%
%   t     current time step (called from ODE solver)
%   z     state                                                    [6 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%   
%   dz    if 3 input arguments:      state derivative wrt time
%
%   Note: It is assumed that the state variables are of the following order:
%         q(1):           [rad]     angle of first joint
%         q(2):           [rad]     angle of second joint
%         q(3):           [rad]     angle of pendulum wrt global z axis
%         dq(1:2):        [rad/s]   angular velocity of joints
%         dq(3):          [rad/s]   angular velocity of pendulum wrt global z axis  
function dz = dynamics_kpdsu(t,z,f1,f2)
%% Code

% set up the system
persistent robot_simul;
if(~isfield(robot_simul, 'A'))
   disp('robot initial construction');
   robot_simul = makeKukaR820_doublep();
end
Vdot_0 = zeros(6,1); Vdot_0(6) = 9.82;
tau = zeros(4,1);
tau(1) = f1(t); 
tau(2) = f2(t);


q = z(1:4);
q(4) = q(4) - q(3);
q(3) = q(3) - q(1) + q(2);

dq = z(5:8);
dq(4) = dq(4) - dq(3);
dq(3) = dq(3) - dq(1) + dq(2);
dz = zeros(8,1);
dz(1:4) = z(5:8);
dz(5:8) = solveForwardDynamics(robot_simul.A,robot_simul.M,q,dq,tau,robot_simul.G, Vdot_0, robot_simul.F);
dz(7)   = dz(5) - dz(6) + dz(7);
dz(8)   = dz(8) + dz(7); 
end