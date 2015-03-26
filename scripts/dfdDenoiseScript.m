%% dfdDenoiseScript

% Download data if needed
%  should we check whether the data has been downloaded already?
% dfdDownloadSampleData

% Load data
for whichSubject = 1:8
    if exist(sprintf(fullfile(dfdRootPath, 'data', 's0%d_sensorData.mat'),whichSubject),'file')
        load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_sensorData.mat'),whichSubject));
    else
        error('Cannot find data file.');
    end
    
    % Preprocess data
    threshold = [0.05 20];
    dataChannels = 1:157;
    sensorData = dfdPreprocessData(sensorData(:,:,dataChannels), threshold);
    
    % Permute sensorData
    sensorData = permute(sensorData, [3,1,2]); % channel x time samples x epoch
    
    % Make design
    if exist(sprintf(fullfile(dfdRootPath, 'data',  's0%d_conditions.mat'),whichSubject),'file')
        load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_conditions.mat'),whichSubject));
    else
        error('Cannot find conditions file.');
    end
    
    %% Make design matrix
    condition_numbers = unique(conditions);
    
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is right (??)
    design(conditions==7,3) = 1; % condition 7 is left (??)    
                                 % condition 3 is blank
    %% Get 'freq' struct to define stimulus locked and broadband frequencies
    %  This struct is needed as input args for getstimlocked and getbroadband
    T = 1;      % epoch length (s)
    fmax = 150; % analyze data up to 150 Hz.
    slF = 12;   % the stimulus flicker frequency (Hz)
    freq = megGetSLandABfrequencies((0:fmax)/T, T, slF/T);
    
    %% Define denoise options for
    % opt.freq            = freq;
    opt.pcstop          = -10;  % denoise with exactly 10 PCs
    opt.preprocessfun   = @hpf; % preprocess data with a high pass filter
    opt.npcs2try        = [];   
    opt.epochGroup      = find(~isnan(squeeze(nanmean(sensorData(:,1,:)))));
    evokedfun           = @(x)getstimlocked(x,freq);
    evalfun             = @(x)getbroadband(x,freq);
    
    
    [results,evalout,denoisedspec,denoisedts] = denoiseData(design,sensorData,evokedfun,evalfun,opt);
        
end
