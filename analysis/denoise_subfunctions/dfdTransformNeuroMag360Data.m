%% Script to transform datasets measured with the NeuroMag360 at CiNet in order to analyze with MEG denoise

%% Add denoise project, meg_utils and fieldtrip
addpath(genpath('/Volumes/server/Projects/MEG/SSMEG/CiNet_SSMEG/code'));
addpath(genpath('~/matlab/git/denoiseproject'));
addpath(genpath('~/matlab/git/meg_utils'))

% ft_path = '/Volumes/server/Projects/MEG/code/fieldtrip';
% meg_add_fieldtrip_paths(ft_path);
%% Define directories
figureDir       = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data'); 

% Define paths and choose your options
whichSubject = 7; % CiNet Subject number
dfdSubject   = 8 + whichSubject; % dfd Subject nr

% Load data
data     = dir(fullfile(dataDir, '*raw*.mat'));
subjData = load(fullfile(dataDir,data(whichSubject).name));

% Cut off the first second for each condition type
for ii = 1:numel(subjData.data_A_all.trial);
    trialGroup = subjData.data_A_all.trial{ii};
    trialGroup = trialGroup(:,1001:end);
    subjData.data_A_all.trial{ii} = trialGroup;
end

for ii = 1:numel(subjData.data_L_all.trial);
    trialGroup = subjData.data_L_all.trial{ii};
    trialGroup = trialGroup(:,1001:end);
    subjData.data_L_all.trial{ii} = trialGroup;
end

for ii = 1:numel(subjData.data_R_all.trial);
    trialGroup = subjData.data_R_all.trial{ii};
    trialGroup = trialGroup(:,1001:end);
    subjData.data_R_all.trial{ii} = trialGroup;
end

% Reshape fieldtrip 'trials' in to epochs 
sensorDataBoth  = trial2epoch(subjData.data_A_all.trial);
sensorDataLeft  = trial2epoch(subjData.data_L_all.trial);
sensorDataRight = trial2epoch(subjData.data_R_all.trial);

% Concatenate sensorData
sensorData = cat(3, sensorDataBoth, sensorDataLeft, sensorDataRight);
sensorData = permute(sensorData, [2 3 1]);

trialNrAll = size(sensorDataBoth,3);

% Make design matrix (check whether all conditions have the same nr of
% epochs)
designB = 3*ones(trialNrAll,1);
designL = 3*ones(trialNrAll+36,1);
designR = 3*ones(trialNrAll+36,1);

for ii = 1:6;
    designB(ii:12:end) = 1;
end

for ii = 1:6;
    designL(ii:12:end) = 7;
end

for ii = 1:6;
    designR(ii:12:end) = 5;
end

conditions = [designB;designL;designR];

% Save dataset & conditions file
save(fullfile(dataDir, sprintf('s%02d_sensorData.mat',dfdSubject+4)),'sensorData');
save(fullfile(dataDir, sprintf('s%02d_conditions.mat',dfdSubject+4)),'conditions');



