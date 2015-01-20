clear all; close all;

%% Script to reproduce Figure 11ABCD S, N, SNR pre-post denoising SL only
% for top ten channels of all subjects
%
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show ...
%
% This script assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

sessionNums        = [1:8];   % You need all subjects for this figure 
                            % Choose a particular number if you would like
                            % a specific subject
plotBb             = false;  % True = Broadband, false = Stim locked                            

inputDataDir       = fullfile(DFDrootpath, 'data'); % Note: New data matrices will 
                                                    % also get stored in the same folder.
whichFun           = 1;     % ToDo: figure out what this means..

conditionNumbers   = 1:6;   % Choose 1:6 to get all conditions: Full, 
                            % left, right (for 'on' 1,2,3 and 'off' 4,5,6)
                            % conditions. 
                            % If you like only a certain conditions e.g.
                            % full on and off define this variable as [1,4]

sensorDataStr      = 'b2';  % Use second version of defining data:
                            % data created by allowing fewer 0 epochs, 
                            % using the new data format (such that design 
                            % matrix and condition order are as they occurred 
                            % in the experiment)
                            
freqDefinition     = 'f';   % f  - uses new definition of ab_i for freq 
                            % (includes more points; excludes 1 pt either side rather than 3)

noisePoolSelection = 'r';   % r - noise pool selection (and pc cutoff, if 
                            % selected by algorithm) by SNR rather than by R2

doHPF              = 'hpf2';% hpf2 - high pass filtered with a sharp cutoff
                            % at 62 Hz and without stimulus harmonics                        
                            
nrPCs              = 'fitfull75'; % fit10 - jump to 10 PCs as cutoff,
                            % In Process: (fitfull10 - don't jump to 10, but denoise all
                            %             channels in between 0 and 10
                            %             PCs.)
                            % fitfull75 - use all channels in noisepool   

xBoots             = 'p1k'; % p1k - bootstrapped 1000x rather than 100x (p100)

fitDataStrBB       = [sensorDataStr freqDefinition noisePoolSelection '_' ...
                        doHPF '_' nrPCs xBoots];
                    
fitDataStrSL       = [sensorDataStr freqDefinition noisePoolSelection 'SL_' ...
                        nrPCs xBoots]; 

saveData         = true;    % Separate matfiles are saved, in order to 
                            % speed up the script if you only want to plot.
saveEpochGroup   = false;   % Epochs can be grouped in a certain order, 
                            % you can save this if you like.
                            
figureDir        = fullfile(DFDrootpath, 'figures');

saveFigures      = true;   % Save figures in the figure folder?

%% Check whether we got our preprocessed data matrices
sessionNums_tmp_BB = [];
sessionNums_tmp_SL = [];

for ii = sessionNums
    % Get session name and top directory
    % This can be modified so that top directory points somewhere else
    dataset = DFDgetdatapaths(ii,conditionNumbers,inputDataDir);
    if ~exist(fullfile(inputDataDir, 'savedProcData', [dataset fitDataStrBB '.mat']),'file');
        sessionNums_tmp_BB(ii) = ii; % Get sessionNumbers for those that dont exist
    end
    if ~exist(fullfile(inputDataDir, 'savedProcData', [dataset fitDataStrSL '.mat']),'file');
        sessionNums_tmp_SL(ii) = ii;
    end
end

if ~isempty(sessionNums_tmp_BB)
    sessionNums_tmp = sessionNums_tmp_BB(sessionNums_tmp_BB ~= 0);
    %% Prepare for denoising (Do we want to separate more part of this preload function?)
    % For example with use of the meg_utils functions?

    % This part just preprocessing and reshaping the input data. 'Bad' channels
    % and 'ok' epochs are defined, and these epochs will get a new value by
    % averaging the surrounding epochs of the same channel. Also, the design
    % matrix will be made, to use later on in the DFDDenoiseAll script.


    % Preload and save data as a 'Session+preprocessing' matlab file. Saved will get a
    % name like "04_SSMEG_04_01_2014b2.mat"

    [sensorData, design, badChannels, conditionNames, okEpochs] = ...
        DFDpreload(sessionNums_tmp, sensorDataStr, saveData, saveEpochGroup, inputDataDir, conditionNumbers);

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
    pcstop10         = false; % To use all PC's in the noise pool (not just number 10)

    resultsBB        = DFDDenoiseWrapper(sessionNums_tmp, [], [], doHpc, [], pcstop10, ...
                                            evalfunToCompute, [], saveDenoiseTs);
end

if ~isempty(sessionNums_tmp_SL)
    sessionNums_tmp = sessionNums_tmp_SL(sessionNums_tmp_SL ~= 0);
    %% Denoise SL for reference
    % We can also can denoise the SL as a check, or use run the function but without 
    % removing any pcs from the data. (In this way you keep the original data).  

    doHpc            = false; % in this case you don't want to high pass the filter including the harmonics
    evalfunToCompute = {'sl'}; % Stimulus locked
    pcstop10         = false; % use all PC's in the noise pool (not just number 10)

    % TODO: With this function you can only denoise with 10 or 75 PCs, if you
    % want to specify any number of PCs yourself, we should change the
    % function. For the SL signal, we will just use the original GLM solution.

    resultsSL        = DFDDenoiseWrapper(sessionNums_tmp, [], [], doHpc, [], pcstop10, ...
                                            evalfunToCompute, [], []);
end
%% Make figure
DFDfiguresnrprepost(sessionNums, conditionNumbers, plotBb, inputDataDir, whichFun, figureDir, saveFigures)