function dfdMakeFigure4()
%% Function to reproduce Figure 4abc (Example data set) for one channel in
% example subject
%
% dfdMakeFigure4()
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
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 



%% Choices to make:
whichSubject    = 1;     % Subject 1 has the example channel.
figureDir       = fullfile(dfdRootPath, 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
                                         
% Define plot colors
colors          = [63, 121, 204; 228, 65, 69; 116,183,74; 127,127,127]/255;
% Define whether to average in log or not
avgLogFlg       = false;

% Load data, design, and get example subject
[data,design,exampleIndex] = prepareData(dataDir,whichSubject,4);

% Define conditions: Full, right, left, off
condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};

condEpochs = {condEpochs1 condEpochs2};

%% Spectrum before and after denoising

fH(1) = plotSpectraPanel4A(data, exampleIndex, condEpochs1, avgLogFlg, colors, saveFigures, figureDir);

fH(2) = plotSpectraPanel4B(data, exampleIndex, condEpochs, avgLogFlg, colors, saveFigures, figureDir);

fH(3) = plotBetaDistributions4C(data, exampleIndex, condEpochs, colors, saveFigures, figureDir); %#ok<NASGU>


