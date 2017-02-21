function dfdMakeFigure5()
%% Function to reproduce Figure 6AB spectra and beta distributionsof full 
% field vs blank, pre-post denoising for one example subject (S1).
%
% dfdMakeFigure5()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show subject's the broadband spectra before and after denoising from 60-150 Hz.
% Bootstrapped fullfield signal (mean across bootstraps) and noise component (std across bootstraps) 
% and the difference between the two distributions are plotted before and
% after denoising. 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:
whichSubject    = 1;     % Subject 1 has the example channel.
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figNum          = 5;     % To call corresponding part of subfunction
                                         
% Define plot colors
colors          = dfdGetColors(4);
% Define whether to average in log or not
avgLogFlg       = false;

% Load data, design, and get example subject
[data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,5);

% Define conditions: Full, right, left, off
condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};

condEpochs = {condEpochs1 condEpochs2};


%% Plot spectrum before and after denoising
fH(1) = plotSpectraPanelBeforeAfterDenoising(data, exampleIndex, exampleChannel, condEpochs, avgLogFlg, colors, saveFigures, figureDir, figNum, 'Figure 6a');

%% Plot signal and noise distributions before and after denoising
fH(2) = plotBetaDistributions(data, exampleIndex, exampleChannel, condEpochs, colors, saveFigures, figureDir, figNum, 'Figure 6b'); %#ok<NASGU>
