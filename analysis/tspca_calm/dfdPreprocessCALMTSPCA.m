function dfdPreprocessCALMTSPCA(whichSubjects)

%% Function to denoise sqd data with CALM or TSPCA algorithms
%
% dfdPreprocessCALMTSPCA(whichSubjects)
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This function will load denoise raw sqd data with CALM (TSPCA with shift
% = 1) or TSPCA denoising algorithms and resaves files. Default parameters
% of TSPCA are based on work of A. de Cheveigne and J. Simon (2007). J.
% Neurosci. Methods.

% "Noise fields measured by reference magnetometers are optimally filtered
% and subtracted from brain channels. The filters (one per reference/brain
% sensor pair) are obtained by delaying the reference signals,
% orthogonalizing them to obtain a basis, projecting the brain sensors onto
% the noise-derived basis, and removing the projections to obtain clean
% data." (DOI: 10.1016/j.jneumeth.2007.06.003).

% For more info about time shifted PCA, see:
% http://www.isr.umd.edu/Labs/CSSL/simonlab/Denoising.html
% http://lumiere.ens.fr/Audition/adc/meg/

% This function depends on meg utils and fieldtrip toolboxes
addpath(genpath('~/matlab/git/meg_utils'));

%% -------------------------------------
% -------------- Load data -------------
% --------------------------------------

dataDir             = '/Volumes/server/Projects/MEG/NYUAD/';
subjectPths         = dir(fullfile(dataDir,'*SSMEG*'));
howToDenoise        = 'CALM'; 
sizeOfBlocks        = 20000;
nearNeighbors       = 0;
excludedChannels    = [21,34];
doWeZeroSaturatedChannels = 'no';
channelForSaturatingChannels = 234; %168;
doWeUseVibrationalChannels = 'no';

switch howToDenoise
    case 'TSPCA'
        postFix = howToDenoise;
        shifts  = [-100:100];
    
    case 'CALM'
        postFix = howToDenoise;
        shifts  = 1;
end


%% -------------------------------------
% ------- STRUCTURE THE RAW DATA -------
% --------------------------------------


%% Load in data

for whichSubject = whichSubjects
  
    dataPth = fullfile(dataDir, subjectPths(whichSubject).name, 'raw'); 
    subjPth    = '*SSMEG*';
    d = dir(fullfile(dataPth, subjPth));
    
    rawFile   = fullfile(dataDir, subjectPths(whichSubject).name, 'raw', d(2).name);
    
    % Get data [timepoints x channels]
    data = sqdread(rawFile);
    
    rawFileName = sprintf('%s_01.con',  d(1).name(1:end-4));
    denoisedFileName = sprintf('%s_%s.con',  d(1).name(1:end-4), postFix);
    
    copyfile(fullfile(dataPth,rawFileName),fullfile(dataPth,denoisedFileName));
    sqdwrite(fullfile(dataPth,rawFileName),fullfile(dataPth,denoisedFileName), 'data', data);

%% -------------------------------------
% ------------ Perform TSPCA -----------
% --------------------------------------

% Time-shift PCA for environmental denoising
sqdDenoise(sizeOfBlocks, shifts, nearNeighbors, rawFile, excludedChannels, ...
            doWeZeroSaturatedChannels, channelForSaturatingChannels, ...
            doWeUseVibrationalChannels, fullfile(dataPth,denoisedFileName));
    
end

