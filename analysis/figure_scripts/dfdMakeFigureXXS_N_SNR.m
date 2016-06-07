function dfdMakeFigure11()
%% Function to reproduce Figure 11ABCD S, N, SNR pre-post denoising SL only
% for top ten channels of all subjects
%
% dfdMakeFigure11()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show the signal component (mean across bootstraps), the
% noise component (std across bootstraps) and signal-to-noise ratio (SNR)
% for stimulus-locked component for three example subjects (Figure 11ABC).
% This figure will show changes in SNR before and after denoising for all
% sessions in the three conditions (full, left, right; Figure 11D).
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function.

%% Choices to make:
whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'exampleAnalysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch');% Where to save images?
saveFigures          = true;     % Save figures in the figure folder?
exampleSessions      = [3,4,5];  % Helena's plot contained subjects [5,6,9]
condColors           = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
dataAll              = [];
figureNumber         = 11;

%% Load data for all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,figureNumber);
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
end

%% S, N, and SNR shown separately, before versus after denoising with 10
%% PCs for stimulus-locked signal. For 3 example sessions - Fig. 11A,B,C
fH = plotSNRPrePostPanelABC(dataAll, exampleSessions, condColors,figureDir,saveFigures,figureNumber); %#ok<NASGU>

%% Plot changes in SNR before and after denoising, showing all sessions
%% together - Fig. 11D
fH = plotSNRPrePostPanelD(dataAll,whichSubjects,condColors,figureDir,saveFigures,figureNumber); %#ok<NASGU>