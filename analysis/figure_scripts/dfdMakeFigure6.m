function dfdMakeFigure6()
%% Function to reproduce Figure 7ABCD S, N, SNR pre-post denoising
% for top ten channels of all subjects
%
% dfdMakeFigure6()
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
whichSubjects    = [1:8];     % Subject 1 has the example channel.
exampleSessions = 1;
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figureNumber    = 6;
                                         
% Define plotting parameters
colors          = dfdGetColors(3);
axmax           = 10;    % How far out do you want to plot the number of PCs

dataAll = [];
for whichSubject = whichSubjects
    fprintf('Load data subject %d \n', whichSubject);
    % Load data, design, and get example subject
    [dataAll{whichSubject},design{whichSubject}] = prepareData(dataDir,whichSubject,7);
end

%% Plot SNR vs number of PCs change for all channels 

fH(1) = plotSNRvsPCsExampleSubjectsPanel6A(dataAll,exampleSessions,colors,axmax,figureDir,saveFigures, 'Figure 6a');

fH(2) = plotSNRvsPCsAllSubjectsPanel7B(dataAll,colors,axmax,figureDir,saveFigures, 'Figure 6b');

fH(3) = plotSNRPrePostPanel(dataAll,whichSubjects,colors,figureDir,saveFigures,figureNumber, 'Figure 6c');

