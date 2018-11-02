%% dynamics_kp_nop.m
% *Summary:* Implements ths ODE for simulating the kuka 6 dof manipulator
%
%    function dz = dynamics_kp_nop(t,z,f)
%
%
% *Input arguments:*
%
%   t     current time step (called from ODE solver)
%   z     state                                                    [4 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%   
%   dz    if 3 input arguments:      state derivative wrt time
%
%   Note: It is assumed that the state variables are of the following order:
%         q(2):           [rad]     angle of joints
%         dq(2):          [rad/s]   angular velocity of joints

function dz = dynamics_NOGP_learning(t,z,f1,f2)
%% Code

% set up the system
persistent robot_simul;
if(~isfield(robot_simul, 'A'))
   disp('robot initial construction');
   robot_simul = makeKukaR820_planar();
end
Vdot_0 = zeros(6,1); Vdot_0(6) = 9.82;
tau = zeros(3,1);
tau(1) = f1(t); 
tau(2) = f2(t);


q = z(1:3);
q(3) = q(3) - q(1) + q(2);
dq = z(4:6);
dq(3) = dq(3) - dq(1) + dq(2);
dz = zeros(6,1);
dz(1:3) = z(4:6);
dz(4:6) = solveForwardDynamics(robot_simul.A,robot_simul.M,q,dq,tau,robot_simul.G, Vdot_0, robot_simul.F);
dz(6)   = dz(4) - dz(5) + dz(6);
end