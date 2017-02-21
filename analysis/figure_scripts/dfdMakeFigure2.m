function dfdMakeFigure2()
%% Function to reproduce Figure 2abc (Example data set) for one channel in
% example subject
%
% dfdMakeFigure2()
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show (a) an the Fourier transformed amplitudes of an example 
% channel for baseline and full flicker condition. (b) High frequency
% portion before and after high pass filtering. (c) Distribution of broadband 
% and stimulus-locked SNR values (Fiu.  
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This function assumes that data is downloaded with the dfdDownloaddata
% function. 


%% Choices to make:
whichSubject    = 1;     % Subject 1 has the example channel.
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figNum          = 2;     % To call corresponding part of subfunction
                                         
% Define plot colors
colors          = dfdGetColors(4);
% Define whether to average in log or not
avgLogFlg       = false;

% Load data, design, and get example subject
[data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,figNum);

% Define conditions: Full, right, left, off
condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};

condEpochs = {condEpochs1 condEpochs2};

%% Spectrum before and after denoising

fH(1) = plotSpectraPanelBeforeDenoising(data, exampleIndex, exampleChannel, condEpochs1, avgLogFlg, colors, saveFigures, figureDir, 'Figure 4a');

fH(2) = plotSpectraPanelBeforeAfterDenoising(data, exampleIndex, exampleChannel, condEpochs, avgLogFlg, colors, saveFigures, figureDir, figNum, 'Figure 4b');

fH(3) = plotBetaDistributions(data, exampleIndex, exampleChannel, condEpochs, colors, saveFigures, figureDir,figNum, 'Figure 4c'); %#ok<NASGU>


