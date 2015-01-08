%% Script to make Figure 1c  
% 
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show .. 
%
%

%% Choices to make:
download_data = true;

% These are the standard options, and can be changed here (This needs to go on top of the script):
sessionNums         = 1:8;
sensorDataStr       = 'b2'; % some sort of data type, we have explain the other options as well
saveData            = true;
saveEpochGroup      = false;
saveData            = true;
savefigures         = false;
conditionNumbers    = 1:6; % Full, left, right (for 'on' and 'off')


%% Get Root paths
project_path = DFDrootpath;

%% Download data

% Usually you want this to be saved in the project_path, but we don't want
% to save large amounts of data in our github repository.

if download_data; 
    save_path = '/Users/winawerlab/Desktop/'; DFDdownloaddata(save_path); 
else
    save_path = '/Users/winawerlab/Desktop/'; 
end

%% Prepare for denoising (Do we want to separate more part of this preload function?)
inputDataDir = save_path;

% Preload and save data
[sensorData, design, badChannels, conditionNames, okEpochs] = ...
    DFDpreload(sessionNums, sensorDataStr, saveData, saveEpochGroup, inputDataDir, conditionNumbers);

%% Denoise BB
dohpc = true;
evalfunToCompute = {'bb'};
results_BB = DFDDenoiseAll(sessionNums, [], [], dohpc, [], [], evalfunToCompute, [], []);

%% Get SL for reference
dohpc = false;
evalfunToCompute = {'sl'};

results_SL = DFDDenoiseAll(sessionNums, [], [], dohpc, [], [], evalfunToCompute, [], []);

%% Plot spatial map

inputDataDir = '~/Desktop/';
fitDataStr = 'b2fr_hpf2_fit10p1k';
%fitDataStr = 'b2frSL_fit10p1k';
whichfun   = 1;
figuredir = 'manuscript_figs/figure_spatialmap';

%% Make spatial map figures:
% 1. example subject spatial map: SL(F,L,R)m Broadband: (F,R,L) = Fig.8
% 2. Plot all subjects - full condition, post minus pre - Fig. 9
% 3. All subjects - right minus left - Fig. 9b (before and after separately)

DFDfigurespatialmap(sessionNums, conditionNumbers,inputDataDir, fitDataStr, whichfun, figuredir,savefigures); % Get the example subject, broadband signal before / after denoising
