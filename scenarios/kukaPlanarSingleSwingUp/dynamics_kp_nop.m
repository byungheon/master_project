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

function dz = dynamics_kp_nop(t,z,f1,f2)
%% Code

% set up the system
persistent robot_kuka_planar_nop;
if isempty(robot_kuka_planar_nop)
   disp('robot initial construction');
   robot_kuka_planar_nop = makeKukaR820_planar_prior;
end
Vdot_0 = zeros(6,1); Vdot_0(6) = 9.82;
tau = zeros(2,1);
tau(1) = f1(t);
tau(2) = f2(t);


q = z(1:2);
dq = z(3:4);
dz = zeros(4,1);
dz(1:2) = dq;
dz(3:4) = solveForwardDynamics(robot_kuka_planar_nop.A,robot_kuka_planar_nop.M,q,dq,tau(:),robot_kuka_planar_nop.G, Vdot_0, robot_kuka_planar_nop.F);

end