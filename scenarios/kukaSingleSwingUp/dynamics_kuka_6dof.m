%% dynamics_kuka_6dof.m
% *Summary:* Implements ths ODE for simulating the kuka 6 dof manipulator
%
%    function dz = dynamics_kssu(t,z,f)
%
%
% *Input arguments:*
%
%   t     current time step (called from ODE solver)
%   z     state                                                    [12 x 1]
%   f     (optional): force f(t)
%
% *Output arguments:*
%   
%   dz    if 3 input arguments:      state derivative wrt time
%
%   Note: It is assumed that the state variables are of the following order:
%         q(6):           [rad]     angle of joints
%         dq(6):          [rad/s]   angular velocity of joints

function dz = dynamics_kuka_6dof(t,z,f1,f2,f3,f4,f5,f6)
%% Code

% set up the system
persistent robot_kuka_6dof;
if(~isfield(robot_kuka_6dof, 'A'))
   disp("robot initial construction");
   robot_kuka_6dof = makeKukaR820_prior;
end
Vdot_0 = zeros(6,1); Vdot_0(6) = 0;
tau = zeros(6,1);
tau(1) = f1(t); 
tau(2) = f2(t);
tau(3) = f3(t);
tau(4) = f4(t);
tau(5) = f5(t);
tau(6) = f6(t);

q = z(1:6);
dq = z(7:12);
dz = zeros(12,1);
dz(1:6) = dq;
dz(7:12) = solveForwardDynamics(robot_kuka_6dof.A,robot_kuka_6dof.M,q,dq,tau(:),robot_kuka_6dof.G, Vdot_0, robot_kuka_6dof.F);

end