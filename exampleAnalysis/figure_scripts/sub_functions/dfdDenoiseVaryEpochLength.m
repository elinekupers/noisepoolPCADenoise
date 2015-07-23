function dfdDenoiseVaryEpochLength(subjects)

% dfdDenoiseVaryEpochLength(subjects)
%
% INPUTS:
% subjects  : Number of subjects one would like to denoise
%
% DESCRIPTION: Function to denoise multiple MEG visual steady
% state data sets, varying in epoch length and number of PCs projected out.
% All results will be saved in one big cell called 'AllResults' and can
% be then be used to reproduce the two supplementary figures from the
% paper the paper:
%
% AUTHORS. YEAR. TITLE. JOURNAL.


% ------------------------------------------------------------------------
% --------------- Define variables depending on how to denoise -----------
% ------------------------------------------------------------------------

% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
dataChannels        = 1:157;
use3Channels        = false;

%% Get frequencies to define stimuluslocked and asynchronous broadband power
% Exclude all frequencies that are close to a multiple of the
% stimulus-locked frequency
f        = 1:150; % Used frequencies
tol      = 1;     % Tolerance for finding stimulus-locked frequency
sl_freq  = 12;    % Stimulus-locked frequency
sl_drop  = f(mod(f, sl_freq) <= tol | mod(f, sl_freq) > sl_freq - tol);
   
% Exclude all frequencies that are close to a multiple of the
% line noise frequency (60 Hz)
ln_drop   = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lf_drop = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);

keep_frequencies    = @(x) x(ab_i);

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,sl_freq); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,keep_frequencies,1000);  % function handle to compuite broadband with a sample rate of 1 kHz

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
        
        opt.epochgroup = epochGroup;
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
    
% ------------------------------------------------------------------------
% -------------------- Denoise and save the data -------------------------
% ------------------------------------------------------------------------
    
    fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data', 's0%d_denoisedData_varyEpochLength_NrPCs'),whichSubject);   
        
    parsave([fname '_bb.mat'], 'allResults', allResults, 'epochDurs', epochDurs, 'npcs', npcs);
    fprintf('data saved:%s\n', fname);
    
end

  
end