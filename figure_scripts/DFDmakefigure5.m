clear all; close all;

%% Script to reproduce Figure 5 (Spatialmap) from example subject 
%
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This script assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

% Denoising options:
% Standard options are applied, and can be changed here.

% TODO: Put the other options also here. Now the standard options are
% defined within the functions itself. 

sessionNums         = 1;        % Choose 1:8 if you would like all the subjects.
                                % But subject 1 is the example subject.
sensorDataStr       = 'b2';     % Some sort of data type. TODO: Explain and 
                                % put other types here as well.
saveData            = true;     % Separate matfiles are saved, in order to 
                                % speed up the script if you only want to plot.
saveEpochGroup      = false;    % Epochs can be grouped in a certain order, 
                                % you can save this if you like.
saveFigures         = true;    % Save figures in the figure folder?
conditionNumbers    = 1:6;      % Choose 1:6 to get all conditions: Full, 
                                % left, right (for 'on' 1,2,3 and 'off' 4,5,6)
                                % conditions. 
                                % If you like only a certain conditions e.g.
                                % full on and off define this variable as [1,4]

fitDataStrBB        = 'b2fr_hpf2_fit10p1k'; % this string is needed to get the correct saved beta values
fitDataStrSL        = 'b2frSL_fit10p1k';
fitDataStr         = 'b2fr_hpf2_fit10p1k';
whichFun            = 1; % I don't exactly know what the difference is betweeen function 1 and 2
figureDir           = fullfile(DFDrootpath, 'figures');

% Get directory where data is saved.
inputDataDir = fullfile(DFDrootpath, 'data'); % Note: New data matrices will 
                                              % also get stored in the same folder.

% We only check for subject 1, since we only need this subject to reproduce plot
[sessionDir] = DFDgetdatapaths(sessionNums,conditionNumbers,inputDataDir);

if ~exist(fullfile(inputDataDir, 'savedProcData', [sessionDir fitDataStrBB '.mat']),'file');
    
    %% Prepare for denoising (Do we want to separate more part of this preload function?)
    % For example with use of the meg_utils functions?

    % This part just preprocessing and reshaping the input data. 'Bad' channels
    % and 'ok' epochs are defined, and these epochs will get a new value by
    % averaging the surrounding epochs of the same channel. Also, the design
    % matrix will be made, to use later on in the DFDDenoiseAll script.
    
    % Preload and save data as a 'Session+preprocessing' matlab file. Saved will get a
    % name like "04_SSMEG_04_01_2014b2.mat"
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

    doHpc            = true; % High pass filter data 
    evalfunToCompute = {'bb'}; % Broadband
    saveDenoiseTs    = true; % You need denoised ts to make the spectrum figure.

    resultsBB        = DFDDenoiseWrapper(sessionNums, [], [], doHpc, [], [], ...
                                            evalfunToCompute, [], saveDenoiseTs);
end

if ~exist(fullfile(inputDataDir, 'savedProcData', [sessionDir fitDataStrSL '.mat']),'file');
    %% Denoise SL for reference
    % We can also can denoise the SL as a check, or use run the function but without 
    % removing any pcs from the data. (In this way you keep the original data).  

    doHpc            = false; % in this case you don't want to high pass the filter including the harmonics
    evalfunToCompute = {'sl'}; % Stimulus locked

    % TODO: With this function you can only denoise with 10 or 75 PCs, if you
    % want to specify any number of PCs yourself, we should change the
    % function. For the SL signal, we will just use the original GLM solution.

    resultsSL        = DFDDenoiseWrapper(sessionNums, [], [], doHpc, [], [], ...
                                            evalfunToCompute, [], []);  
end
                                                            

%% Make spatial map figures:
% The DFDfigurespatialmap function can make three types of spatial maps.
% Define with the first term which figure you want to make: 
% Option 5: Reproduce Figure 5. Example subject spatial map: SL(F,L,R) & Broadband: (F,R,L)
% Option 8: Reproduce Figure 8. Plot full condition, post minus pre
%                               denoising
% Option 9: Reproduce Figure 9. Plot right minus left stimulation, SL and BB 
%                               pre and post denosing separately)

% TODO: add plotting functions to the aux folder, so it will not depend on Fieldtrip toolbox on server 
DFDfigurespatialmap(5, sessionNums, conditionNumbers,inputDataDir, fitDataStrBB, whichFun, figureDir,saveFigures);



