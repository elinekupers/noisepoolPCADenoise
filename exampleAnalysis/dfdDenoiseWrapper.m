function dfdDenoiseWrapper(subjects, howToDenoise)
%
% dfdDenoiseWrapper(subjects, howToDenoise)
%
% INPUTS:
% subjects  : Number of subjects one would like to denoise
% howToDenoise: 1, 2, or 3, meaning:
%                   1: denoise with exactly 10 PC regressors
%                   2: denoise with each of 0 to 10 PC regressors
%                   3: various controls
%               [default=1]
%
% DESCRIPTION: Wrapper function to denoise multiple MEG visual steady
% state data sets. The results of the denoising are written to file and can
% be then be used to reproduce several of the published figures from the
% paper the paper:
%
% AUTHORS. YEAR. TITLE. JOURNAL.

% ------------------------------------------------------------------------
% --------------------------- Check options ------------------------------
% ------------------------------------------------------------------------

% Check for how to denoise the data
if notDefined('howToDenoise'), howToDenoise = 1; end

% Check for data, download data if needed
if isempty(dir(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 's0*.mat')));
    error('No data were found. Use dfdDownloadSampleData')
end

% ------------------------------------------------------------------------
% --------------- Define variables depending on how to denoise -----------
% ------------------------------------------------------------------------

% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
dataChannels        = 1:157;
use3Channels        = false;

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
freq                = megGetSLandABfrequencies(0:150, 12, 60, 1);

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,freq);  % function handle to compuite broadband

% Define options for denoising that are equal for each type of denoising
opt.resampling      = {'boot','boot'};
opt.pcselmethod     = 'snr';

switch howToDenoise % Define denoise other parameters (see denoisedata.m)
    case 1 % Denoise with exactly 10 PC regressors
        opt.pcchoose          = -10;   % Take 10 PCs
        opt.npcs2try          = [];    
        optsl                 = opt;
        optbb                 = opt;
        optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
        nrControlModes        = 0;
        postFix               = '';
        
    case 2 % Denoise with each of 0 to 10 PC regressors
        opt.pcchoose          = 1.05;   % Get threshold for optimal nr of PCs
        opt.npcs2try          = 10;     % loop through 10
        optsl                 = opt;
        optbb                 = opt;
        optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
        nrControlModes        = 0;
        postFix               = 'full';
        
    case 3 % Denoise with various control modes
        opt.pcchoose          = 10;     % Take 10 PCs
        opt.npcs2try          = [];     
        optsl                 = opt;
        optbb                 = opt;
        optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
        nrControlModes        = 1:5;   % All control modes
        postFix               = 'control';       
end

% ------------------------------------------------------------------------
% ------------------ Load and denoise data per subject -------------------
% ------------------------------------------------------------------------

for whichSubject = subjects
    % ------------------ Load data and design ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    % ------------------ Make design matrix ------------------------------
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is left field
    design(conditions==7,3) = 1; % condition 7 is right field
    % condition 3 is blank
    
    % ------------------ Preprocess data ---------------------------------
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
    
    % ---- Remove bad channels and bad epochs from data and conditions ---
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design     = design(~badEpochs,:);
    
    % ------------------ Permute sensorData for denoising ----------------
    sensorData = permute(sensorData, [3 1 2]);  
    
% ------------------------------------------------------------------------
% -------------------- Denoise and save the data -------------------------
% ------------------------------------------------------------------------
    
    % ------------------ Denoise for broadband analysis ------------------
    for nrControl = nrControlModes;                                                                                               
        if (0 <= nrControl) && (nrControl <= 4)
            optbb.pccontrolmode   = nrControl;
            postFix               = sprintf('control%d',nrControl);
            [results,evalout]     = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
        elseif nrControl == 5
            optbb.pccontrolmode   = 0;
            optbb.npoolmethod     = {'r2','n',size(sensorData,1)};
            postFix               = sprintf('control%d',nrControl);
            [results,evalout]     = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
        end
        
        % ------------------ Define file name ---------------------------
        if use3Channels
            fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data', ['s0%d_denoisedData_' postFix '_w3chan']),whichSubject);
        else
            fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data', ['s0%d_denoisedData_' postFix]), whichSubject);   
        end
        
        % ----------------- Save denoised broadband data -----------------
        parsave([fname '_bb.mat'], 'results', results, 'evalout', evalout, 'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)
        
        % ----------- Denoise and save stimulus-locked analysis ----------
        [results,evalout] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);
        parsave([fname '_sl.mat'], 'results', results, 'evalout', evalout, 'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optsl)
    end
end
