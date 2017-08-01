function dfdMakeFigureS1()
%% Function to reproduce supplementary Figure 3 SNR pre-post denoising
% when removing pcs from a synthetic data set that contains 5 or 10 bases
% functions
%
% dfdMakeFigureS1()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show the simulated broadband response before and after denoising.
% Bootstrapped fullfield signal (mean across bootstraps) and noise component (std across bootstraps)
% and the difference between the two distributions are plotted before and
% after denoising.
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function.

%% Choices to make:
whichSubject    = 99;     % Subject 99 has the synthetic dataset.
exampleSessions = 1;      % this session contains 10 basis functions
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figureNumber    = 'SF1';
ylimsSF3A       = [-5, 80];
ylimsSF3B       = [-1, 28];

% Define plotting parameters
colors          = dfdGetColors(3);
axmax           = 75;    % How far out do you want to plot the number of PCs

fprintf('Load data subject %d \n', whichSubject);
% Load data, design, and get example subject
dataAll = prepareData(dataDir,whichSubject,figureNumber);

%% Create power spectrum, before denoising


%% Broadband before and after denoising


%% Plot SNR vs number of PCs change for all channels
fH(2) = plotSNRvsPCsAllSubjectsPanel7B({dataAll},colors,axmax,figureDir,saveFigures,ylimsSF3B, 'Figure SF1D');

fH(3) = plotNoisepool(99);

fH(4) = dfdMakeFigure9(whichSubject); for ax = 2:4:16; set(fH(3).Children(ax),'CLim', [-20 20], 'YTick', [-20 -10 0 10 20]); end

