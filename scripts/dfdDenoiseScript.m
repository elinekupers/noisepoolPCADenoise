%% dfdDenoiseScript

% Add paths
dfdAddpaths

% Download data if needed
dfdDownloadSampleData

% Load data
for whichSubject = 1:8
    if exist(sprintf(fullfile(dfdRootPath, 'data', 'raw', 's0%d_sensorData.mat'),whichSubject),'file')
        load(sprintf(fullfile(dfdRootPath, 'data','raw', 's0%d_sensorData.mat'),whichSubject));
    else
        disp('Error. Cannot find data file.');
    end
    
    % Preprocess data
    threshold = [0.05 20];
    dataChannels = 1:157;
    sensorData = dfdPreprocessData(sensorData(:,:,dataChannels), threshold);
    
    % Permute sensorData
    sensorData = permute(sensorData, [3,1,2]); % channel x time samples x epoch
    
    % Make design
    if exist(sprintf(fullfile(dfdRootPath, 'data', 'raw', 's0%d_conditions.mat'),whichSubject),'file')
        load(sprintf(fullfile(dfdRootPath, 'data','raw', 's0%d_conditions.mat'),whichSubject));
    else
        disp('Error. Cannot find conditions file.');
    end
    
    %% Make design matrix
    condition_numbers = unique(conditions);
    
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1;
    design(conditions==5,2) = 1;
    design(conditions==7,3) = 1;
    
    %% Get index frequency range to define SL and BB
    T = 1; fmax = 150; slF = 12;
    freq = megGetSLandABfrequencies((0:fmax)/T, T, slF/T);
    
    %% Define denoise options for
    opt.freq            = freq;
    opt.pcstop          = -10;
    opt.preprocessfun   = @hpf;
    opt.npcs2try        = 10;
    opt.epochGroup      = ~isnan(squeeze(sensorData(1,1,:)));
    
    evokedfun           = @(x)getstimlocked(x,freq);
    evalfun             = @(x)getbroadband(x,freq);
    
    
    [results,evalout,denoisedspec,denoisedts] = denoiseData(design,sensorData,evokedfun,evalfun,opt);
        
end

% Plot data
dfdMakeAllFigures;
