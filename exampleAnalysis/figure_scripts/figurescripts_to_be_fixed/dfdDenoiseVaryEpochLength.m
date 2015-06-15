function dfdDenoiseVaryEpochLength(subjects)

% denoise as a function of denoising epoch duration 



% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
dataChannels        = 1:157;
use3Channels        = false;

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
freq = megGetSLandABfrequencies(0:150, 1, 12);

% Define options for denoising
opt.resampling        = {'boot','boot'};
opt.pcselmethod       = 'snr';
opt.preprocessfun     = @hpf;  % preprocess data with a high pass filter for broadband analysis
opt.npoolmethod       = {'snr','n',75};

% Define functions to define noise pool and signal of interest
evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband

% Get different epoch lengths and npcs to denoise with
epochDurs             = [1,3,6,12,24,36,72,inf];
npcs                  = [5,10:10:70];

allResults = [];

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
    design = design(~badEpochs,:);
    
    % ------------------ Permute sensorData for denoising ----------------
    sensorData = permute(sensorData, [3 1 2]);  

    % ------------------ Loop over Epoch length & nr PCs -----------------
    for ii = 1:length(epochDurs) % iterate through epoch duration fo doing PCA
        
        if isinf(epochDurs(ii))
            epochGroup = ones(size(sensorData,1),1);
        else
            epochGroup = megEpochGroup(~badEpochs,epochDurs(ii),0);
        end
        
        opt.epochGroup = epochGroup;
        fprintf('epochDur = %d; %d epochs\n', epochDurs(ii), max(epochGroup));
        
        for jj = 1:length(npcs) % iterate through number of PCs projected out 
            opt.pcstop = -npcs(jj);
            fprintf('\tnpcs = %d\n', npcs(jj));
            
            if jj == 1
                [results] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
                noisepooldef = results.noisepool;
            else
                [results] = denoisedata(design,sensorData,noisepooldef,evalfun,opt);
            end
            allResults{ii,jj} = results;
            
            clear results;
        end
    end
    
    fname = fullfile(sprintf(dfdRootPath,'exampleAnalysis','data', 's0%d_denoisedData_varyEpochLength_NrPCs',whichSubject));   
        
    parsave([fname '_bb.mat'], 'allResults', allResults, 'epochDurs', epochDurs, 'npcs', npcs);
    fprintf('data saved:%s\n', fname);
    
end

  
end