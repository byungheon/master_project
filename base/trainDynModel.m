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

in_list = [];
for i = 1:size(targets_temp,1)-1
   tmpsign = sign(targets_temp(i,plant.jointi)).*sign((inputs_temp(i,plant.jointi+length(plant.jointi))+inputs_temp(i+1,plant.jointi+length(plant.jointi)))/2);
   tmpmag  = abs(targets_temp(i,plant.jointi)/dt)./ abs((inputs_temp(i,plant.jointi+length(plant.jointi))+inputs_temp(i+1,plant.jointi+length(plant.jointi)))/2);
   tmp = tmpsign.*(tmpmag <2.0).*(tmpmag>0.5);
   if  tmp>0
%        disp('----------------------');
%        disp(targets_temp(i,plant.jointi)/dt);
%        disp((inputs_temp(i,plant.jointi+length(plant.jointi))+inputs_temp(i+1,plant.jointi+length(plant.jointi)))/2);
       in_list = [in_list i];
   end
end
if length(in_list) > 100
   disp('Random Sampling from obtained dynamic data');
   new_list = randsample(length(in_list),100);
   in_list = in_list(new_list);
end
x(end - initial_length +1:end,:) = [];
y(end - initial_length +1:end,:) = [];
xx = xx(in_list,:);
yy = yy(in_list,:);
x = [x;xx];
y = [y;yy];
targets_temp    = targets_temp(in_list,:);
inputs_temp     = inputs_temp(in_list,:);

disp(['Obtained Dynamics Data: ' num2str(length(in_list)) '/' num2str(initial_length)]);
if isfield(dynmodel, 'inputs')
    dynmodel.inputs     = [dynmodel.inputs;inputs_temp];
    dynmodel.targets    = [dynmodel.targets;targets_temp]; 
else
    dynmodel.inputs     = inputs_temp;
    dynmodel.targets    = targets_temp;
end


if (isfield(dynmodel,'model') && ~strcmp(dynmodel.model,'PILCO'))
    delta           = zeros(size(dynmodel.targets));
    delta(:,plant.jointi)    = dynmodel.inputs(:,plant.jointi+length(plant.jointi)) * dt;
    for i = 1:size(dynmodel.inputs,1)
        delta(i,plant.jointi+length(plant.jointi)) = dt * solveForwardDynamics(dynmodel.robot.A,dynmodel.robot.M,dynmodel.inputs(i,plant.jointi)',dynmodel.inputs(i,plant.jointi+length(plant.jointi))',dynmodel.inputs(i,end-Du+1:end)',dynmodel.robot.G, dynmodel.Vdot0, dynmodel.robot.F);
    end
    dynmodel.targets = dynmodel.targets - delta;
end

dynmodel = dynmodel.train(dynmodel, plant, trainOpt);  %  train dynamics GP

% display some hyperparameters
Xh = dynmodel.hyp;     
% noise standard deviations
disp(['Learned noise std: ' num2str(exp(Xh(end,:)))]);
% signal-to-noise ratios (values > 500 can cause numerical problems)
disp(['SNRs             : ' num2str(exp(Xh(end-1,:)-Xh(end,:)))]);