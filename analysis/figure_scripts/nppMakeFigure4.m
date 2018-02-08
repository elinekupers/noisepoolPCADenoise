function nppMakeFigure4()
%% Function to reproduce Figure 4abc (Example data set) for one channel in
% example subject
%
% nppMakeFigure4()
% 
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will show (a) an the Fourier transformed amplitudes of an
% example channel for baseline and full flicker condition. (b) High
% frequency portion before and after high pass filtering. (c) Distribution
% of broadband and stimulus-locked SNR values (signal-blue and noise-gray
% (upper) and difference-black (lower)).
%
% This function assumes that data is downloaded with the nppDownloaddata
% function.


%% Choices to make:
whichSubject    = 1;     % Subject 1 has the example channel.
figureDir       = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figNum          = 4;     % To call corresponding part of subfunction
                                         
% Define plot colors
colors          = nppGetColors(4);
% Define whether to average in log or not
avgLogFlg       = false;

% Load data, design, and get example subject
[data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,figNum);

% Define conditions: Full, right, left, off
condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};

condEpochs = {condEpochs1 condEpochs2};

%% Spectrum before and after denoising

fH(1) = plotSpectraPanelBeforeDenoising(data, exampleIndex, exampleChannel, condEpochs1, avgLogFlg, colors, saveFigures, figureDir, 'Figure 4A');

fH(2) = plotSpectraPanelBeforeAfterDenoising(data, exampleIndex, exampleChannel, condEpochs, avgLogFlg, colors, saveFigures, figureDir, figNum, 'Figure 4B');

fH(3) = plotBetaDistributions(data, exampleIndex, exampleChannel, condEpochs, colors, saveFigures, figureDir,figNum, 'Figure 4C'); %#ok<NASGU>


