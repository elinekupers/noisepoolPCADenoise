function dfdMakeFigureSXX2GridNPCvsNoisepool()

%% Function to reproduce an obsolete supplementary figure 

% dfdMakeFigureSXX2GridNPCvsNoisepool()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show a grid for the three conditions (full, left, right), 
% where the difference in SNR before and after denoising is plotted for number
% of PCs removed versus number of channels in the noisepool. 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function and then analyzed by dfdDenoiseNPCvsNoisepool(1:8).


%% Choices to make:
whichSubjects        = [1:8]; 
dataDir              = fullfile(dfdRootPath, 'analysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'analysis', 'figures');% Where to save images?
saveFigures          = true;     % Save figures in the figure folder?
dataAll              = [];
figureNumber         = 'SFXX2';
npools               = [5,10:10:140];
npcs                 = [5,10:10:130];

%% Load data for all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,figureNumber);
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
end

%% Plot difference in SNR (post-pre) as a function of number of channels in
%  noise pool and number of PCs removed
fH = plotNPCvsNoisepoolSNR(whichSubjects,dataAll,npools,npcs,saveFigures,figureDir);
