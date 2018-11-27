
%% cost graph
% (optional) Plot predicted immediate costs (as a function of the time steps
if exist('plotting', 'var') && isfield(plotting, 'verbosity') ...
    && plotting.verbosity > 0
  if ~ishandle(3); figure(3); else set(0,'CurrentFigure',3); end
  clf(3); errorbar(0:dt:(H*dt),fantasy.mean{j},2*fantasy.std{j}); drawnow;
  xlabel('Time[s]','FontSize', 24); ylabel('Immediate cost','FontSize', 24);
end
hold on; plot(0:dt:((length(realCost{J+j})-1)*dt),realCost{J+j},'r','LineWidth',2); drawnow;
lgd = legend('Predicted','Real');
lgd.FontSize = 30;
%% sequential cost graph
number_pilco = 10;
number_mine = 6;
pilco_cost = zeros(1,number_pilco);
mine_cost = zeros(1,number_mine);
for index = 1:number_pilco
    name = 'C:\Users\BH\Desktop\졸논프로젝트\pilco_dynamics\data\Success\PILCO\pilco_5g_11085_12trials\KukaSwingUp_5g_11085_3313_PILCO_' + string(index) + '_H200.mat';
%     name = 'C:\Users\BH\Desktop\졸논프로젝트\pilco_dynamics\data\Success\PILCO\pilco_5g_11085_12trials\KukaSwingUp_5g_11085_3313_PILCO_' + string(index) + '_H200.mat';
    load(name);
    pilco_cost(index) = sum(realCost{J+j});
end

for index = 1:number_mine
    name = 'C:\Users\BH\Desktop\졸논프로젝트\pilco_dynamics\data\Success\MINE\mine_5g_11085_6trials\KukaSwingUp_5g_11085_3313_MINE_' + string(index) + '_H200.mat';
%     name = 'C:\Users\BH\Desktop\졸논프로젝트\pilco_dynamics\data\Success\MINE\12pmipmi_4_0024_6trials_8_8\KukaSwingUp_12pmipmi_4_0024_088MINE_' + string(index) + '_H200.mat';

    load(name);
    mine_cost(index) = sum(realCost{J+j});
end

figure(4);
plot(1:number_pilco,pilco_cost,'b-o');
hold on;
plot(1:number_mine,mine_cost,'r-o');
xlabel('Rollouts[#]','FontSize', 24); ylabel('Long-term cost','FontSize', 24);
lgd = legend('PILCO','Ours');
lgd.FontSize = 15;
%% theta prediction graph
% close all;
load('F:\bh\data\KukaSwingUp_12pmipmi_4_0024_MINE_8_H200.mat');
% load('C:\Users\BH\Desktop\pilco_dynamics\data\Success\MINE\12pmipmi_4_0024_8trials_25_8\KukaSwingUp_12pmipmi_4_002_4_MINE_8_H200')

q1_pred = M{j}(1,1:end-1);
q2_pred = M{j}(2,1:end-1);
q3_pred = M{j}(6,1:end-1);

q1 = xx(:,1)';
q2 = xx(:,2)';
q3 = xx(:,6)';

figure(10);
plot(0:dt:199*dt, q1_pred,'b','LineWidth',1.5);
hold on;
plot(0:dt:199*dt, q1,'r','LineWidth',1.5);
xlabel('Time[s]','FontSize', 18); ylabel('Angle of the first joint [rad]','FontSize', 18);
lgd = legend('predicted','real','Location','southeast');
lgd.FontSize = 30;

figure(11);
plot(0:dt:199*dt, q2_pred,'b','LineWidth',1.5);
hold on;
plot(0:dt:199*dt, q2,'r','LineWidth',1.5);
xlabel('Time[s]','FontSize', 18); ylabel('Angle of the second joint [rad]','FontSize', 18);
lgd = legend('predicted','real','Location','southeast');
lgd.FontSize = 30;

figure(12);
plot(0:dt:199*dt, q3_pred,'b','LineWidth',1.5);
hold on;
plot(0:dt:199*dt, q3,'r','LineWidth',1.5);
xlabel('Time[s]','FontSize', 18); ylabel('Angle of the pendulum [rad]','FontSize', 18);
lgd = legend('predicted','real','Location','southeast');
lgd.FontSize = 30;

%%
clear all;
clc;
close all;

load('C:\Users\BH\Desktop\졸논프로젝트\pilco_dynamics\data\Success\MINE\12pmipmi_4_0024_6trials_8_8\KukaSwingUp_12pmipmi_4_0024_088MINE_6_H200.mat');
T = 10;
HH = ceil(T/dt);           % prediction steps (optimization horizon)
[xx, yy, realCost{j+J}, latent{j} ss] = ...
  rollout(gaussian(mu0, S0), policy, HH, plant, cost, dynmodel);
% disp(xx);                           % display states of observed trajectory
if ~ishandle(3); figure(3); else set(0,'CurrentFigure',3); end
hold on; plot(1:length(realCost{J+j}),realCost{J+j},'r'); drawnow;

%%
q_draw = latent{j}(:,[plant.jointi plant.angi(3)]);
q_draw(:,3) = q_draw(:,3) - q_draw(:,1) + q_draw(:,2);

appVisualizeKUKA_video(q_draw', dt);