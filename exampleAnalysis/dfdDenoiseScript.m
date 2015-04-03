%% dfdDenoiseScript

% Description: Script to denoise 8 example MEG visual steady state data
% sets. The results of the denoising are used to make several figures for
% the paper:
%
% AUTHORS. YEAR. TITLE. JOURNAL.

opt.verbose = true;

% Check for data, download data if needed
if isempty(fullfile(dfdRootPath, 'data'));
    error('No data were found. Use dfdDownloadSampleData')
end

% Load and denoise data, one subject at a time
for whichSubject = 1:8
    % ******* Load data and design *****************
    load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_sensorData.mat'),whichSubject));
    load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_conditions.mat'),whichSubject));
    
    % ******* Preprocess data **********************
    varThreshold        = [0.05 20];
    badChannelThreshold = 0.2;
    badEpochThreshold   = 0.2;
    dataChannels        = 1:157;
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, opt.verbose);
    
    % ******** Make design matrix *****************
    condition_numbers = unique(conditions);
    
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is right (??)
    design(conditions==7,3) = 1; % condition 7 is left (??)
    % condition 3 is blank
    
    % Remove bad channels and bad epochs from data and conditions
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design = design(~badEpochs,:);
    
    % Get 'freq' struct to define stimulus locked and broadband frequencies
    %  This struct is needed as input args for getstimlocked and getbroadband
    T = 1;      % epoch length (s)
    fmax = 150; % analyze data up to 150 Hz.
    slF = 12;   % the stimulus flicker frequency (Hz)
    freq = megGetSLandABfrequencies((0:fmax)/T, T, slF/T);
        
    % ********* Denoise the data ********************

    % Permute sensorData for denoising
    sensorData = permute(sensorData, [3 1 2]);
    
    %   Denoise for broadband analysis
    opt.pcchoose        = -10;   % denoise with exactly 10 PCs
    opt.preprocessfun   = @hpf;  % preprocess data with a high pass filter
    evokedfun           = @(x)getstimlocked(x,freq);
    evalfun             = @(x)getbroadband(x,freq);
    
    [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,opt); %#ok<ASGLU>
    save(sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData_bb.mat'),whichSubject), ...
        'results','evalout','denoisedspec','denoisedts','badChannels','badEpochs','-v7.3');
    
    %   Denoise for stimulus-locked analysis
    opt.preprocessfun   = [];
    evalfun             = @(x)getstimlocked(x,freq);
    
    [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    save(sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData_sl.mat'),whichSubject), ...
        'results','evalout','denoisedspec','denoisedts','badChannels','badEpochs','-v7.3');
    
end
