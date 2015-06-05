function dfdMakeFigure6()
%% Function to reproduce Figure 6A SNR against PCs removed for 3 example subjects 
% and Figure 6b, SNR for top ten channels against PCs removed for all
% subjects.
%
% dfdMakeFigure6()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show 2 ways of plotting SNR values against the number of PCs
% removed from the data. In Figure 6A, SNR values of all channels from the 
% three example subjects (corresponding to 6,7,8) will be plotted against
% the number of PCs removed. The median is plotted in the color of the
% condition. In Figure 6B, only the median of all channels for SNR against
% number of PCs removed, for all subjects, for the three conditions. This
% function needs all 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:
whichSubjects       = 1:8;             
figureDir           = fullfile(dfdRootPath, 'exampleAnalysis', 'figures'); % Where to save images?
dataDir             = fullfile(dfdRootPath, 'exampleAnalysis', 'data');    % Where to save data?
saveFigures         = true;         % Save figures in the figure folder?
exampleSessions     = [3,4,5];
condColors          = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
axmax               = 10;           % how far to go out on the x-axis

%% Load data from all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,6);    
    dataAll{whichSubject} = {data,design,exampleIndex};   %#ok<AGROW>
end

%% SNR increase as a function of number of PCs removed, 3 example sessions - Fig. 6A
fH = plotSNRvsPCsExampleSubjectsPanel6A(dataAll,exampleSessions,condColors,axmax,figureDir,saveFigures); %#ok<NASGU>

%% SNR increase as a function of number of PCs removed, all sessions - Fig. 6B
fH = plotSNRvsPCsAllSubjectsPanel6B(dataAll,condColors,axmax,figureDir,saveFigures); %#ok<NASGU>

