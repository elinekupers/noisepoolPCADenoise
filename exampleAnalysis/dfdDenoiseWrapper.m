function dfdDenoiseWrapper(subjects)

% Description: Wrapper function to denoise multiple MEG visual steady
% state data sets. The results of the denoising are written to file and can
% be then be used to reproduce several of the published figures from the
% paper the paper:
%
% AUTHORS. YEAR. TITLE. JOURNAL.

% Check for data, download data if needed
if isempty(fullfile(dfdRootPath, 'data'));
    error('No data were found. Use dfdDownloadSampleData')
end

% preprocessing parameters (see dfdPreprocessData)
varThreshold        = [0.05 20];  
badChannelThreshold = 0.2;       
badEpochThreshold   = 0.2;
dataChannels        = 1:157;
use3Channels        = true;

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
freq = megGetSLandABfrequencies(0:150, 1, 12);

% denoise parameters (see denoisedata.m)
optsl.pcchoose        = -10;   % denoise with exactly 10 PCs for stim`ulus locked
optbb.pcchoose        = -10;   % denoise with exactly 10 PCs for broadband
optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband

% Load and denoise data, one subject at a time
parfor whichSubject = subjects
    % ******* Load data and design *****************
    tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    % ******** Make design matrix *****************
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is right (??)
    design(conditions==7,3) = 1; % condition 7 is left (??)
    % condition 3 is blank
    
    % ******* Preprocess data **********************
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
    
    % Remove bad channels and bad epochs from data and conditions
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design = design(~badEpochs,:);

    % Permute sensorData for denoising
    sensorData = permute(sensorData, [3 1 2]);
    
    % ********* Denoise the data ********************    
    
    %   Denoise for broadband analysis    
    [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
    
    if use3Channels
        fname = sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData_w3chan'),whichSubject);
    else
        fname = sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData'),whichSubject);
    end
    
    parsave([fname '_bb.mat'], 'results', results, 'evalout', evalout, ...
        'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
        'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)        
    
    %   Denoise for stimulus-locked analysis
    [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);
    parsave([fname '_sl.mat'], 'results', results, 'evalout', evalout, ...
        'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
        'badChannels', badChannels, 'badEpochs', badEpochs,  'opt', optsl)        
    
end
