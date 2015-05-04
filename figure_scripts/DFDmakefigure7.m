function dfdMakeFigure7()

%% Script to reproduce Figure 7ABCD S, N, SNR pre-post denoising
% for top ten channels of all subjects
% subjects.
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show ...
%
% This script assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

whichSubjects        = 1:8; 
dataDir              = fullfile(dfdRootPath, 'data');
figureDir            = fullfile(dfdRootPath, 'figures');
saveFigures          = true;   % Save figures in the figure folder?
exampleSessions      = [3,4,5];  % Helena's plot contained subjects [5,6,9]
condColors           = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
dataAll              = [];
%% Load data

for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);

    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,7);
    
    dataAll{whichSubject} = {data,design,exampleIndex};
    
    
end
%% S, N, and SNR shown separately, before versus after denoising with 10
%% PCs. For 3 example sessions - Fig. 7A,B,C

fH = plotSNRPrePostPanel7ABC(dataAll, exampleSessions, condColors,figureDir,saveFigures)


%% Plot changes in SNR before and after denoising, showing all sessions
%% together - Fig. 7D
fH = plotSNRPrePostPanel7D(dataAll,whichSubjects,condColors,figureDir,saveFigures)

