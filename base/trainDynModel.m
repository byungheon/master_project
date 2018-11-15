%% trainDynModel.m
% *Summary:* Script to learn the dynamics model
%
% Copyright (C) 2008-2013 by 
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modification: 2013-05-20
%
%% High-Level Steps
% # Extract states and controls from x-matrix
% # Define the training inputs and targets of the GP
% # Train the GP 

%% Code

% 1. Train GP dynamics model
initial_length = size(xx,1);
Du = length(policy.maxU); Da = length(plant.angi); % no. of ctrl and angles
xaug = [xx(:,dyno) xx(:,end-Du-2*Da+1:end-Du)];     % x augmented with angles
inputs_temp = [xaug(:,dyni) xx(:,end-Du+1:end)];     % use dyni and ctrl
targets_temp = yy(:,dyno);
targets_temp(:,difi) = targets_temp(:,difi) - xx(:,dyno(difi));

% marginalize gp dynamics
joint_temp = ss(:,[jointi length(jointi) + jointi]);

if (isfield(dynmodel,'model') && ~strcmp(dynmodel.model,'PILCO'))
    local.dynamics    = @dynamics_kp_nop;
    local.OPTIONS     = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
    local.ctrlfcn     = str2func('zoh_local');   
    local.par.dt      = dt; local.par.delay = 0; local.par.tau = dt;
    local.u0          = cell(1,length(jointi));
    local.delta       = zeros(1,length(difi));
    for i = 1:size(targets_temp,1)
        local.tau = inputs_temp(i,end-Du+1:end)';
        
%         [~, local.y] = ode45(@(t,input)dynamics_kp_nop_not(t,input,local.tau(1),local.tau(2)), [0 dt/2 dt], joint_temp(i,:)', local.OPTIONS);
%         
%         local.qdotdelt  = local.y(3,jointi) - joint_temp(i,jointi);
%         local.qddotdelt = local.y(3,jointi + length(jointi)) - joint_temp(i,jointi + length(jointi));
%         
%         local.delta(plant.jointi)           = local.qdotdelt;
%         local.delta(jointi+length(jointi))  = local.qddotdelt;
        
        local.qddot   = solveForwardDynamics(dynmodel.robot.A,dynmodel.robot.M,joint_temp(i,jointi)',joint_temp(i,jointi + length(jointi))',local.tau,dynmodel.robot.G,dynmodel.Vdot0, dynmodel.robot.F);
        local.delta(plant.jointi)           = (joint_temp(i,jointi + length(jointi)) + local.qddot' * dt/2) * dt;
        local.delta(jointi+length(jointi))  = local.qddot' * dt;

        targets_temp(i,:)                   = targets_temp(i,:) - dynmodel.dynratio * local.delta;
    end
    clear local;
end

% random sampling from gathered dynamics data
if size(targets_temp,1) > 100
   disp('Random Sampling from obtained dynamic data');
   in_list = randsample(size(targets_temp,1),100);
end
x(end - initial_length +1:end,:) = [];
y(end - initial_length +1:end,:) = [];
xx = xx(in_list,:);
yy = yy(in_list,:);
x = [x;xx];
y = [y;yy];
targets_temp    = targets_temp(in_list,:);
inputs_temp     = inputs_temp(in_list,:);

% augment current dynamics data to previous data
disp(['Obtained Dynamics Data: ' num2str(length(in_list)) '/' num2str(initial_length)]);
if isfield(dynmodel, 'inputs')
    dynmodel.inputs     = [dynmodel.inputs;inputs_temp];
    dynmodel.targets    = [dynmodel.targets;targets_temp]; 
else
    dynmodel.inputs     = inputs_temp;
    dynmodel.targets    = targets_temp;
end

% train dynamics data
dynmodel = dynmodel.train(dynmodel, plant, trainOpt);  %  train dynamics GP

% display some hyperparameters
Xh = dynmodel.hyp;     
% noise standard deviations
disp(['Learned noise std: ' num2str(exp(Xh(end,:)))]);
% signal-to-noise ratios (values > 500 can cause numerical problems)
disp(['SNRs             : ' num2str(exp(Xh(end-1,:)-Xh(end,:)))]);

function u = zoh_local(f, t, par) % **************************** zero-order hold
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