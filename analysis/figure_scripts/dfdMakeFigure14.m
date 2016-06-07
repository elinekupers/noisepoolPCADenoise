function dfdMakeFigure14
%% Function to reproduce Figure 13 (Spatialmap) across CiNet dataset subjects
%
% dfdMakeFigure5AcrossSubjects(whichSubjects,figureDir,dataDir,saveFigures,threshold)
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

%% Get denoising results for raw (no preprocessing) data
whichSubjectsRaw    = 1:8;                % Raw
dfdMakeFigure14AcrossSubjects(whichSubjectsRaw,figureDir,dataDir,saveFigures,threshold)

%% Get denoising results for CALM preprocessed data
whichSubjectsCALM   = 21:28;              % CALM
dfdMakeFigure14AcrossSubjects(whichSubjectsCALM,figureDir,dataDir,saveFigures,threshold)

%% Get denoising results for TSPCA preprocessed data
whichSubjectsTSPCA  = 29:36;              % TSPCA
dfdMakeFigure14AcrossSubjects(whichSubjectsTSPCA,figureDir,dataDir,saveFigures,threshold)
