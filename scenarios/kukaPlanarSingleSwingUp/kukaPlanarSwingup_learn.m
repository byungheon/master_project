%% kukaPlanarSingleSwingup_learn.m
% *Summary:* Script to learn a controller for the cart-doube-pendulum
% swingup
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-27
%
%% High-Level Steps
% # Load parameters
% # Create J initial trajectories by applying random controls
% # Controlled learning (train dynamics model, policy learning, policy
% application)

%% Code

% 1. Initialization
clear all; close all;
settings_kpssu;                 % load scenario-specific settings
basename = 'KukaPlanarSingleSwingUp_'; % filename used for saving data
% 2. Initial J random rollouts
for jj = 1:J
  [xx, yy, realCost{jj}, latent{jj}] = ...
    rollout(gaussian(mu0, S0), struct('maxU',policy.maxU), H, plant, cost, dynmodel);
  x = [x; xx]; y = [y; yy];       % augment training sets for dynamics model
  if plotting.verbosity > 0;      % visualization of trajectory
    draw_rollout_kpssu; 
  end
end
% writetable(array2table(x), 'x.txt','Delimiter',',');
% writetable(array2table(y), 'y.txt','Delimiter',',');
% writetable(array2table(realCost{jj}), 'realCost.txt','Delimiter',',');
% writetable(array2table(latent{jj}), 'latent.txt','Delimiter',',');
% x = dlmread('x.txt');
% y = dlmread('y.txt');
% xx= x;
% yy = y;
% realCost{1} = dlmread('realCost.txt');
% latent{1} = dlmread('latent.txt');
% jj = 1;

mu0Sim(odei,:) = mu0; S0Sim(odei,odei) = S0;
mu0Sim = mu0Sim(dyno); S0Sim = S0Sim(dyno,dyno);

% 3. Controlled learning (N iterations)
for j = 1:N
  tic
  trainDynModel;   % train (GP) dynamics model
  learnPolicy;     % learn policy
  applyController; % apply controller to system
  disp(['controlled trial # ' num2str(j)]);
  if plotting.verbosity > 0      % visualization of trajectory
    draw_rollout_kpssu; 
  end 
  toc
end
