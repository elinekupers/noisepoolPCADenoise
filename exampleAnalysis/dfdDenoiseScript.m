%% dfdDenoiseScript

% Description

saveResults = true;

% Check for data, download data if needed
if isempty(fullfile(dfdRootPath, 'data'));
    error('No data were found. Use dfdDownloadSampleData')
end

% Load data
for whichSubject = 1%:8
    % Load data and design
    load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_sensorData.mat'),whichSubject));
    load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_conditions.mat'),whichSubject));
    
    % Preprocess data
    varThreshold        = [0.05 20];
    badChannelThreshold = 0.2;
    badEpochThreshold   = 0.2;
    dataChannels        = 1:157;
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold);              
            
    %% Make design matrix
    condition_numbers = unique(conditions);
    
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is right (??)
    design(conditions==7,3) = 1; % condition 7 is left (??)    
                                 % condition 3 is blank
    
    % Remove bad channels and bad epochs from data and conditions
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design = design(~badEpochs,:);
    
    %% Get 'freq' struct to define stimulus locked and broadband frequencies
    %  This struct is needed as input args for getstimlocked and getbroadband
    T = 1;      % epoch length (s)
    fmax = 150; % analyze data up to 150 Hz.
    slF = 12;   % the stimulus flicker frequency (Hz)
    freq = megGetSLandABfrequencies((0:fmax)/T, T, slF/T);
    
    %% Permute sensorData for denoising
    sensorData = permute(sensorData, [3 1 2]);
    
    %% Define denoise options for
    opt.pcstop          = -10;  % denoise with exactly 10 PCs
    opt.preprocessfun   = @hpf; % preprocess data with a high pass filter
    opt.npcs2try        = [];   
    opt.badepochs       = badEpochs;    
    opt.badchannels     = badChannels;
    evokedfun           = @(x)getstimlocked(x,freq);
    evalfun             = @(x)getbroadband(x,freq);
      
    [bbresults,bbevalout,bbdenoisedspec,bbdenoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
  
    % Do it for SL as well
    evokedfun           = @(x)getstimlocked(x,freq);
    evalfun             = @(x)getstimlocked(x,freq);

    [slresults,slevalout,sldenoisedspec,sldenoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    
    if saveResults
        save(sprintf(fullfile(dfdRootPath,'data','s0%d_preprocData.mat'),whichSubject), ...
                'sensorData');
        save(sprintf(fullfile(dfdRootPath,'data','s0%d_bbdenoisedData.mat'),whichSubject), ...
                'bbresults','bbevalout','bbdenoisedspec','bbdenoisedts','badChannels','badEpochs');
        save(sprintf(fullfile(dfdRootPath,'data','s0%d_sldenoisedData.mat'),whichSubject), ...
                'slresults','slevalout','sldenoisedspec','sldenoisedts','badChannels','badEpochs');
        save(sprintf(fullfile(dfdRootPath,'data','s0%d_preprocDesign.mat'),whichSubject),'design');
    end
        
end
