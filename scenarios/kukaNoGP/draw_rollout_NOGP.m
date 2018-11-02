%% draw_rollout_kssu.m
% *Summary:* Script to draw a trajectory of the kuka-single-swingup system and
%
% Copyright (C) 2008-2013 by
% Marc Deisenroth, Andrew McHutchon, Joe Hall, and Carl Edward Rasmussen.
%
% Last modified: 2013-03-27
%
%% High-Level Steps
% # For each time step, plot the observed trajectory and the predicted
% means and covariances of the Cartesian coordinates of the tips of both
% pendulums

%% Code

if exist('j','var') && ~isempty(M{j})
  q_draw = latent{j}(:,[plant.jointi]);
  text1 = ['trial # ' num2str(j+J) ', T=' num2str(H*dt) ' sec'];
  text2 = ['total experience (after this trial): ' num2str(H*dt*j) ' sec'];
else
  q_draw = latent{jj}(:,[plant.jointi]);
  text1 = ['(random) trial # ' num2str(jj) ', T=' num2str(H*dt) ' sec'];
  text2 = ['total experience (after this trial): ' num2str(H*dt*jj) ' sec'];
end
q_draw(:,3) = q_draw(:,3) - q_draw(:,1) + q_draw(:,2);
appVisualizeKUKA_planar(q_draw', cost, text1, text2, dt);
