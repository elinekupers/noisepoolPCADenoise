%% Script to reproduce Figure 6A SNR against PCs removed for 3 example subjects 
% and Figure 6b, SNR for top ten channels against PCs removed for all
% subjects.
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show ... 
%
% This script assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

sessionNums      = 1:8;     % You need all subjects for this figure 
                            % Choose a particular number if you would like
                            % a specific subject
plotBb           = true;    % True = Broadband, false = Stim locked                            
inputDataDir     = fullfile(DFDrootpath, 'data'); % Note: New data matrices will 
                                                  % also get stored in the same folder.
whichFun         = 1;       % ToDo: figure out what this means..
conditionNumbers = 1:6;     % Choose 1:6 to get all conditions: Full, 
                            % left, right (for 'on' 1,2,3 and 'off' 4,5,6)
                            % conditions. 
                            % If you like only a certain conditions e.g.
                            % full on and off define this variable as [1,4]
sensorDataStr    = 'b2';
                                              
fitDataStrBB     = 'b2fr_hpf2_fitfull75p1k'; % this string is needed to get the correct saved beta values
fitDataStrSL     = 'b2frSL_fitfull75p1k'; % this string is needed to get the correct saved beta values
saveData         = true;    % Separate matfiles are saved, in order to 
                            % speed up the script if you only want to plot.
saveEpochGroup   = false;   % Epochs can be grouped in a certain order, 
                            % you can save this if you like.
figureDir        = fullfile(DFDrootpath, 'figures');
saveFigures      = false;   % Save figures in the figure folder?


for ii = sessionNums
    % Get session name and top directory
    % This can be modified so that top directory points somewhere else
    dataset = DFDgetdatapaths(ii,conditionNumbers,inputDataDir);
    if ~exist(fullfile(inputDataDir, 'savedProcData', [dataset fitDataStrBB '.mat']),'file');
        preloaddataBB = true;
    end
    if ~exist(fullfile(inputDataDir, 'savedProcData', [dataset fitDataStrSL '.mat']),'file');
        preloaddataSL = true;
    end
    if preloaddataBB && preloaddataSL == true
        break
    end
end


if preloaddataBB
    %% Prepare for denoising (Do we want to separate more part of this preload function?)
    % For example with use of the meg_utils functions?

    % This part just preprocessing and reshaping the input data. 'Bad' channels
    % and 'ok' epochs are defined, and these epochs will get a new value by
    % averaging the surrounding epochs of the same channel. Also, the design
    % matrix will be made, to use later on in the DFDDenoiseAll script.


    % Preload and save data as a 'Session+preprocessing' matlab file. Saved will get a
    % name like "04_SSMEG_04_01_2014b2.mat"

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

    doHpc            = true; % High pass filter data 
    evalfunToCompute = {'bb'}; % Broadband
    saveDenoiseTs    = true; % You need denoised ts to make the spectrum figure.

    resultsBB        = DFDDenoiseWrapper(sessionNums, [], [], doHpc, [], [], ...
                                            evalfunToCompute, [], saveDenoiseTs);
end

if preloaddataSL
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
%% Make figure
DFDfiguresnrversuspcs(sessionNums, plotBb, inputDataDir, whichFun, figureDir, saveFigures)