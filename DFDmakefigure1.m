%% Script to reproduce Figure 5 (Spatialmap) from...
% 
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show .. 
%
%

%% Choices to make:

% Denoising options:
% Standard options are applied, and can be changed here.

% TODO: Put the other options also here. Now the standard options are
% defined within the functions itself. 

sessionNums         = 1;        % Choose 1:8 if you would like all the subjects 
sensorDataStr       = 'b2';     % Some sort of data type. TODO: Explain and put other types here as well.
saveData            = true;     % Separate matfiles are saved, in order to speed up the script if you only want to plot.
saveEpochGroup      = false;    % Epochs can be grouped in a certain order, you can save this if you like.
savefigures         = false;    % Save figures in the subject folder?
conditionNumbers    = 1:6;      % Choose 1:6 to get all conditions: Full, left, right (for 'on' 1,2,3 and 'off' 4,5,6) conditions. 
%                                   If you like only a certain conditions e.g.
%                                   full on and off define this variable as [1,4]


%% Get Root paths
project_path = DFDrootpath; % folder with all the project functions

%% Download data

% Usually you want this to be saved in the project_path, but we don't want
% to save large amounts of data in our github repository.

if download_data; 
    save_path = '/Users/winawerlab/Desktop/'; DFDdownloaddata(save_path); 
else
    save_path = '/Users/winawerlab/Desktop/'; 
end

%% Prepare for denoising (Do we want to separate more part of this preload function?)
% For example with use of the meg_utils functions?

% This part just preprocessing and reshaping the input data. 'Bad' channels
% and 'ok' epochs are defined, and these epochs will get a new value by
% averaging the surrounding epochs of the same channel. Also, the design
% matrix will be made, to use later on in the DFDDenoiseAll script.

inputDataDir = save_path; % all our data is stored on the Desktop for now.

% Preload and save data as a 'Session+preprocessing type' 
% (Example:04_SSMEG_04_01_2014b2.mat)
[sensorData, design, badChannels, conditionNames, okEpochs] = ...
    DFDpreload(sessionNums, sensorDataStr, saveData, saveEpochGroup, inputDataDir, conditionNumbers);

%% Denoise BB

% In the DFDDenoiseAll function, we will conduct another preprocessing
% step, before we actually denoise the data. This preprocessing step is schematically 
% shown in Figure 2a,b,c of the manuscript. First, define SL and BB frequencies.
% Second, calculate the values of these two signal types, for every channel and epoch. 
% 
% The next analysis steps are described in Figure 3a,b,c and are part of
% the denoisedata function:
% First, define a noise pool with the channels having the lowest SNR values
% for the SL signal.
% Second, one can choose to high pas filter the broadband signals, this is especially
% needed for denoising the BB signals. For the SL, you don't want to do
% this, because the 12 Hz SL and its harmonics will be filtered out. (Which
% is something you actually want to keep).
% Conduct PCA on every channel and epoch of the noise pool.
% Remove the pre-selected number of PCs from the BB data.

% Denoised data and bootstrapped beta solutions will be saved in a new
% mat-files with the corresponding options in the name.

dohpc = true;
evalfunToCompute = {'bb'}; % Broadband
saveDenoiseTs = true; % You need denoised ts to make the spectrum figure.
results_BB = DFDDenoiseAll(sessionNums, [], [], dohpc, [], [], evalfunToCompute, [], true);

%% Denoise SL for reference
% We can also can denoise the SL as a check, or use run the function but without 
% removing any pcs from the data. (In this way you keep the original data).  


dohpc = false; % in this case you don't want to high pass the filter including the harmonics
evalfunToCompute = {'sl'}; % Stimulus locked

% TODO: With this function you can only denoise with 10 or 75 PCs, if you
% want to specify any number of PCs yourself, we should change the
% function. For the SL signal, we will just use the original GLM solution.

results_SL = DFDDenoiseAll(sessionNums, [], [], dohpc, [], [], evalfunToCompute, [], []);

%% Plot spatial map

% Predefine some variables (TODO: Bring this to the top)
inputDataDir = '~/Desktop/';
fitDataStr = 'b2fr_hpf2_fit10p1k'; % this string is needed to get the correct saved beta values
%fitDataStr = 'b2frSL_fit10p1k';
whichfun   = 1; % I don't exactly know what the difference is betweeen function 1 and 2
figuredir = 'manuscript_figs/figure_spatialmap';

%% Make spatial map figures:
% 1. example subject spatial map: SL(F,L,R)m Broadband: (F,R,L) = Fig.8
% 2. Plot all subjects - full condition, post minus pre - Fig. 9
% 3. All subjects - right minus left - Fig. 10 (?) (before and after separately)

% TODO: add plotting functions to the aux folder, so it will not depend on Fieldtrip toolbox on server 

DFDfigurespatialmap(sessionNums, conditionNumbers,inputDataDir, fitDataStr, whichfun, figuredir,savefigures);

%% Make spectrum figures: 
% These are the same figures as 4a,b,c and make use of data from channel 42.

DFDfigurespectrum(sessionNums, conditionNumbers, inputDataDir, fitDataStr, figuredir, savefigures)