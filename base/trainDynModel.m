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
Du = length(policy.maxU); Da = length(plant.angi); % no. of ctrl and angles
xaug = [x(:,dyno) x(:,end-Du-2*Da+1:end-Du)];     % x augmented with angles
dynmodel.inputs = [xaug(:,dyni) x(:,end-Du+1:end)];     % use dyni and ctrl
dynmodel.targets = y(:,dyno);
dynmodel.targets(:,difi) = dynmodel.targets(:,difi) - x(:,dyno(difi));

in_list = [];
for i = 1:size(dynmodel.targets,1)-1
   tmpsign = sign(dynmodel.targets(i,plant.jointi)).*sign((dynmodel.inputs(i,plant.jointi+length(plant.jointi))+dynmodel.inputs(i+1,plant.jointi+length(plant.jointi)))/2);
   tmpmag  = abs(dynmodel.targets(i,plant.jointi)/dt)./ abs((dynmodel.inputs(i,plant.jointi+length(plant.jointi))+dynmodel.inputs(i+1,plant.jointi+length(plant.jointi)))/2);
   tmp = tmpsign.*(tmpmag <2.0).*(tmpmag>0.5);
   if  tmp>0
       disp('----------------------');
       disp(dynmodel.targets(i,plant.jointi)/dt);
       disp((dynmodel.inputs(i,plant.jointi+length(plant.jointi))+dynmodel.inputs(i+1,plant.jointi+length(plant.jointi)))/2);
       in_list = [in_list i];
   end
end
disp(['Obtained Dynamics Data: ' num2str(length(in_list)) '/' num2str(size(dynmodel.targets,1))]);
dynmodel.inputs     = dynmodel.inputs(in_list,:);
dynmodel.targets    = dynmodel.targets(in_list,:);
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