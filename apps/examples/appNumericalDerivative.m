
%% Inverse 1st Derivative
clear all
close all
clc
robot = makeKukaR820_planar_prior;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
Friction = robot.F;
Vdot0 = [0,0,0,0,0,9.8]';
dt  = 0.001;
x_start = randn();
x   = x_start:dt:(x_start+20*dt);
q       = -2 + 4 * rand(2,1);
qdot    = (-2 + 4 * rand(2,1)) * 2;
qddot   = (-2 + 4 * rand(2,1)) * 5;

tau         = zeros(2,length(x));
dtaudq      = zeros(2,2,length(x));
dtaudqdot   = zeros(2,2,length(x));
dtaudqddot  = zeros(2,2,length(x));

for tau_i = 1:2
    for joint_i = 1:2
        for i =1:length(x)
            q_tmp       = q;
        %     q_tmp(joint_i)    = x(i);
            qdot_tmp    = qdot;
        %     qdot_tmp(joint_i)    = x(i);
            qddot_tmp   = qddot;
            qddot_tmp(joint_i)    = x(i);

            [tau_tmp, T, V, Vdot, F] = solveInverseDynamics(A,M,q_tmp,qdot_tmp,qddot_tmp,G, Vdot0, Friction);
            tau(:,i)  = tau_tmp;
            [dtaudq(:,:,i), dtaudqdot(:,:,i),  dtaudqddot(:,:,i)] = solveInverseDynamicsDerivatives_pilco(A,M,q_tmp,qdot_tmp,G,T,V,Vdot,F,Vdot0,Friction);
        end

        numerical = gradient(tau(tau_i,:),dt);
        numerical = numerical(2:end-1);
        analytic  = dtaudqddot(tau_i, joint_i,:);
        analytic  = analytic(:)';
        analytic = analytic(2:end-1);analytic - numerical
        disp('--------------------------');
    end
end
        

%% Inverse 2nd Derivative
clear all
close all
clc
robot = makeKukaR820_planar_prior;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
Friction = robot.F;
Vdot0 = [0,0,0,0,0,9.8]';
dt  = 0.001;
x_start = randn();
x   = x_start:dt:(x_start+20*dt);
q       = -2 + 4 * rand(2,1);
qdot    = (-2 + 4 * rand(2,1)) * 2;
qddot   = (-2 + 4 * rand(2,1)) * 5;

tau                 = zeros(2,length(x));
dtaudq              = zeros(2,2,length(x));
dtaudqdot           = zeros(2,2,length(x));
dtaudqddot          = zeros(2,2,length(x));
dtaudqdq            = zeros(2,4,length(x));
dtaudqdqdot         = zeros(2,4,length(x));
dtaudqdqddot        = zeros(2,4,length(x));
dtaudqdotdqdot      = zeros(2,4,length(x));
dtaudqdotdqddot     = zeros(2,4,length(x));
dtaudqddotdqddot    = zeros(2,4,length(x));

list_cell = [1,1,1,2,2,3];
cell_dtau = {dtaudq, dtaudqdot, dtaudqddot};
cell_ddtau = {dtaudqdq, dtaudqdqdot, dtaudqdqddot, dtaudqdotdqdot, dtaudqdotdqddot, dtaudqddotdqddot};
for cell_i = 1:6
    for tau_i = 1:2
        for joint_i = 1:2
            for i =1:length(x)
                q_tmp       = q;
                qdot_tmp    = qdot;
                qddot_tmp   = qddot;
                switch cell_i
                    case {1}
                        q_tmp(joint_i)    = x(i);
                    case {2,4}
                        qdot_tmp(joint_i)    = x(i);
                    case {3,5,6}
                        qddot_tmp(joint_i)    = x(i);
                end
                [tau_tmp, T, V, Vdot, F] = solveInverseDynamics(A,M,q_tmp,qdot_tmp,qddot_tmp,G, Vdot0, Friction);
                tau(:,i)  = tau_tmp;
                [cell_dtau{1}(:,:,i), cell_dtau{2}(:,:,i),  cell_dtau{3}(:,:,i)] = solveInverseDynamicsDerivatives_pilco(A,M,q_tmp,qdot_tmp,G,T,V,Vdot,F,Vdot0,Friction);
                [~, ~, ~, cell_ddtau{1}(:,:,i), cell_ddtau{2}(:,:,i), cell_ddtau{3}(:,:,i), cell_ddtau{4}(:,:,i), cell_ddtau{5}(:,:,i), cell_ddtau{6}(:,:,i)] = ...
                    solveInverseDynamicsSecondDerivatives_pilco(A,M,q_tmp,qdot_tmp,G,T,V,Vdot,F,Vdot0,Friction);
            end
            numerical = zeros(2,length(x));
            for i = 1:2
                gradient_tmp = cell_dtau{list_cell(cell_i)}(tau_i,i,:);
                gradient_tmp = gradient_tmp(:)';
                numerical(i,:) = gradient(gradient_tmp,dt);
            end
            numerical = numerical(:,2:end-1);
            analytic  = cell_ddtau{cell_i}(tau_i, [joint_i,joint_i+2],:);
            analytic  = reshape(analytic,[2,size(analytic,3)]);
            analytic = analytic(:,2:end-1);analytic - numerical
            disp('--------------------------');

        end
    end
end

%% Forward 1st Derivative
clear all
close all
clc
robot = makeKukaR820_planar_prior;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
F = robot.F;
Vdot0 = [0,0,0,0,0,9.8]';
dt  = 0.001;
x_start = randn();
x   = x_start:dt:(x_start+20*dt);
q       = -2 + 4 * rand(2,1);
qdot    = (-2 + 4 * rand(2,1)) * 2;
tau     = (-2 + 4 * rand(2,1)) * 10;
qddot   = zeros(2,length(x));
dqddotdq = zeros(2,2,length(x));
dqddotdqdot = zeros(2,2,length(x));
dqddotdtau = zeros(2,2,length(x));


for qddot_i = 1:2
    for input_i = 1:2
        for i =1:length(x)
            q_tmp       = q;
        %     q_tmp(input_i)    = x(i);
            qdot_tmp    = qdot;
        %     qdot_tmp(input_i)    = x(i);
            tau_tmp     = tau;
            tau_tmp(input_i)    = x(i);

            qddot_tmp   = solveForwardDynamics(A,M,q_tmp,qdot_tmp,tau_tmp,G,Vdot0,F);
        %     [tau_inverse, ~, ~, ~, ~] = solveInverseDynamics(A,M,q_tmp,qdot_tmp,qddot_tmp,G, Vdot0, F);
        %     tau_tmp' - tau_inverse'
            qddot(:,i)  = qddot_tmp;
            [dqddotdq(:,:,i), dqddotdqdot(:,:,i), dqddotdtau(:,:,i)] = solveForwardDynamicsDerivatives_pilco(A,M,q_tmp,qdot_tmp,qddot_tmp,G,Vdot0, F);
        end                                                            


        numerical = gradient(qddot(qddot_i,:),dt);
        numerical = numerical(2:end-1);
        analytic  = dqddotdtau(qddot_i, input_i,:);
        analytic  = analytic(:)';
        analytic  = analytic(2:end-1);analytic - numerical
        disp('----------------------------');
    end
end




%% Forward 2nd Derivative
clear all
close all
clc
robot = makeKukaR820_planar_prior;
% robot = makeKukaR820_planar;

n = robot.dof;
A = robot.A;
M = robot.M;
G = robot.G;
Friction = robot.F;
Vdot0 = [0,0,0,0,0,9.8]';
dt = 0.00001;
x_start = randn();
x       = x_start:dt:(x_start+20*dt);
q       = -2 + 4 * rand(n,1);
qdot    = (-2 + 4 * rand(n,1)) * 2;
tau     = (-2 + 4 * rand(n,1)) * 10;
% q = [  -1.0063;
%    -0.7184];
% qdot = [2.6401;
%    -3.6490];
% tau = [-3.5037;
%    -6.5485];

qddot                   = zeros(n,length(x));
dqddotdq                = zeros(n,n,length(x));
dqddotdqdot             = zeros(n,n,length(x));
dqddotdtau              = zeros(n,n,length(x));
dqddotdqdq              = zeros(n,n*n,length(x));
dqddotdqdqdot           = zeros(n,n*n,length(x));
dqddotdqdtau            = zeros(n,n*n,length(x));
dqddotdqdotdqdot        = zeros(n,n*n,length(x));
dqddotdqdotdtau         = zeros(n,n*n,length(x));
dqddotdtaudtau          = zeros(n,n*n,length(x));

list_cell = [1,1,1,2,2,3];
cell_dqddot = {dqddotdq, dqddotdqdot, dqddotdtau};
cell_ddqddot= {dqddotdqdq, dqddotdqdqdot, dqddotdqdtau, dqddotdqdotdqdot, dqddotdqdotdtau, dqddotdtaudtau};
for cell_i = 1:1
for qddot_i = 1:n
    for input_i = 1:n
       for i =1:length(x)
            q_tmp       = q;
            qdot_tmp    = qdot;
            tau_tmp     = tau;
            switch cell_i
                case {1}
                    q_tmp(input_i)    = x(i);
                case {2,4}
                    qdot_tmp(input_i)    = x(i);
                case {3,5,6}
                    tau_tmp(input_i)    = x(i);
            end
            qddot_tmp   = solveForwardDynamics(A,M,q_tmp,qdot_tmp,tau_tmp,G,Vdot0,Friction);
            qddot(:,i)  = qddot_tmp;
            [cell_dqddot{1}(:,:,i), cell_dqddot{2}(:,:,i),  cell_dqddot{3}(:,:,i)] = solveForwardDynamicsDerivatives_pilco(A,M,q_tmp,qdot_tmp,qddot_tmp,G,Vdot0, Friction);
            [cell_dqddot{1}(:,:,i), cell_dqddot{2}(:,:,i),  cell_dqddot{3}(:,:,i), cell_ddqddot{1}(:,:,i), cell_ddqddot{2}(:,:,i), cell_ddqddot{3}(:,:,i), cell_ddqddot{4}(:,:,i), cell_ddqddot{5}(:,:,i), cell_ddqddot{6}(:,:,i)] = ...
                solveForwardDynamicsSecondDerivatives_pilco(A,M,q_tmp,qdot_tmp,qddot_tmp,G,Vdot0,Friction);
       end
        numerical = zeros(n,length(x));

        for i = 1:n
            gradient_tmp = cell_dqddot{list_cell(cell_i)}(qddot_i,i,:);
            gradient_tmp = gradient_tmp(:)';
            numerical(i,:) = gradient(gradient_tmp,dt);
        end
        numerical = numerical(:,2:end-1);
        list_input = [];
        for k = 1:n
            list_input = [list_input input_i+n*(k-1)]; 
        end
        analytic  = cell_ddqddot{cell_i}(qddot_i, list_input,:);
        analytic  = reshape(analytic,[n,size(analytic,3)]);
        analytic  = analytic(:,2:end-1);analytic - numerical 
        disp('--------------------------');

    end
end
end
        