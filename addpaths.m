addpath(genpath('libraries/'));
addpath(genpath('apps/'));
addpath(genpath('models/'));
addpath(genpath('scenarios'));
addpath('base','control','doc','gp','loss','test','util','practice_note');
if exist('data','dir') == 7
   mkdir('data');
end
addpath('data');