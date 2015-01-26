clear all; close all;

%% Script to reproduce Figure 4abc (Example data set) for one channel in example subject
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show (a) an the Fourier transformed amplitudes of an example 
% channel for baseline and full flicker condition. (b) High frequency
% portion before and after high pass filtering. (c) Distribution of broadband 
% and stimulus-locked SNR values.  
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This script assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

sessionNums         = 1;        % Subject 1 has the example channel.
                                % Choose 1:8 if you would like all the subjects.
conditionNumbers    = 1:6;      % Choose 1:6 to get all conditions: Full, 
                                % left, right (for 'on' 1,2,3 and 'off' 4,5,6)
                                % conditions. 
                                % If you like only a certain conditions e.g.
                                % full on and off define this variable as [1,4]
sensorDataStr       = 'b2';     % Use second version of defining data:
                                % data created by allowing fewer 0 epochs, 
                                % using the new data format (such that design 
                                % matrix and condition order are as they occurred 
                                % in the experiment)
                            
freqDefinition      = 'f';      % f  - uses new definition of ab_i for freq 
                                % (includes more points; excludes 1 pt either side rather than 3)

noisePoolSelection  = 'r';      % r - noise pool selection (and pc cutoff, if 
                                % selected by algorithm) by SNR rather than by R2

doHPF               = 'hpf2';   % hpf2 - high pass filtered with a sharp cutoff
                                % at 62 Hz and without stimulus harmonics                        
                            
nrPCs               = 'fit10';  % fit10 - always use 10 PCs as cutoff,
                               % fitfull75 - use all channels in noisepool

xBoots              = 'p1k';     % p1k - bootstrapped 1000x rather than 100x (p100)
                          
fitDataStrBB          = [sensorDataStr freqDefinition noisePoolSelection '_' ...
                        doHPF '_' nrPCs xBoots]; % 'b2fr_hpf2_fit10p1k', this string is needed to get the correct saved beta values
                    
fitDataStrSL          = [sensorDataStr freqDefinition noisePoolSelection 'SL_' ...
                        nrPCs xBoots];  % ''b2frSL_fit10p1k'', this string is needed to get the correct saved beta values
                    
                    
saveData            = true;     % Separate matfiles are saved, in order to 
                                % speed up the script if you only want to plot.
saveEpochGroup      = false;    % Epochs can be grouped in a certain order, 
                                % you can save this if you like.
                                                     
figureDir           = fullfile(DFDrootpath, 'figures');

saveFigures         = true;    % Save figures in the figure folder?

inputDataDir = fullfile(DFDrootpath, 'data'); % Note: New data matrices will 
                                              % also get stored in the same folder.
                                              
% Get session dir for one subject
[sessionDir] = DFDgetdatapaths(sessionNums,conditionNumbers,inputDataDir);

% Check whether the saved data mat file already exists.. Otherwise, make it 
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

    evalfunToCompute = {'bb'}; % Broadband
    saveDenoiseTs    = true; % You need denoised ts to make the spectrum figure.
    doHpc            = true;
    

    resultsBB        = DFDDenoiseWrapper(sessionNums, [], [], doHpc, [], [], ...
                                            evalfunToCompute, [], saveDenoiseTs);
end

if ~exist(fullfile(inputDataDir, 'savedProcData', [sessionDir fitDataStrSL '.mat']),'file');
    %% Denoise SL for reference
    % We can also can denoise the SL as a check, or use run the function but without 
    % removing any pcs from the data. (In this way you keep the original data).  

    evalfunToCompute = {'sl'}; % Stimulus locked

    % TODO: With this function you can only denoise with 10 or 75 PCs, if you
    % want to specify any number of PCs yourself, we should change the
    % function. For the SL signal, we will just use the original GLM solution.

    resultsSL        = DFDDenoiseWrapper(sessionNums, [], [], doHpc, [], [], ...
                                            evalfunToCompute, [], []);  
end


DFDfigurespectrum(sessionNums, conditionNumbers, inputDataDir, fitDataStrBB, figureDir, saveFigures)
