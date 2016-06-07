function dfdMakeFigure13

%% Function to reproduce Figure 13 (Spatialmap) across all subjects of CiNet dataset
%
% dfdMakeFigure13()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the stimulus locked signal, broadband signals before
% using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately).
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function.

%% Choices to make:
figureDir           = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir             = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     	= true;     % Save figures in the figure folder?
threshold           = 0;

%% Get denoising results for raw data
whichSubjectsRaw    = 9:12;                % Raw
dfdMakeFigure13AcrossSubjects(whichSubjectsRaw,figureDir,dataDir,saveFigures,threshold)


%% Get denoising results for tsss data
whichSubjectsTSSS   = [14,16,18,20];       % TSSS
dfdMakeFigure13AcrossSubjects(whichSubjectsTSSS,figureDir,dataDir,saveFigures,threshold)
