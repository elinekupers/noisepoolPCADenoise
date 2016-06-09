function dfdMakeFigureVaryEpochLength()

%% Function to reproduce Supplementary Figure 1 with SNR pre-post difference
% for different epoch lengths and different epoch length against used npcs.
% for top ten channels of all subjects
%
% dfdFigureVaryEpochLength()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show ...
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function and then analyzed by dfdDenoiseVaryEpochLength.


%% Choices to make:
whichSubjects        = 1%[1:8];
dataDir              = fullfile(dfdRootPath, 'analysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'analysis', 'figures');% Where to save images?
saveFigures          = true;     % Save figures in the figure folder?
condColors           = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
dataAll              = [];
figureNumber         = 'SF1';
epochDurs            = [1,3,6,12,24,36,72,1080];
npcs                 = [5,10:10:70];

%% Load data for all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,figureNumber);
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
end

%% Plot difference in SNR (post-pre) as a function of denoising epoch duration
fH = plotEpochLengthVersusSNR(whichSubjects,dataAll,epochDurs,condColors,saveFigures,figureDir);

%% Plot surface of difference in SNR (post-pre) as a function of epoch duration and
%% number of PCs removed. 
% Some specs:
%  - Noisepool selection by SNR,
%  - Highpass filtered,
%  - 10 PCs removed
%  - bootstrapped 100 x. 

fH = plotEpochLengthVersusNPCsVersusSNR(whichSubjects,dataAll,epochDurs,npcs,saveFigures,figureDir); %#ok<*NASGU>
