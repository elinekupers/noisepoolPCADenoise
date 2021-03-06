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
removeFirstEpoch    = true;

%% Get frequencies to define stimuluslocked and asynchronous broadband power
% Exclude all frequencies that are close to a multiple of the
% stimulus-locked frequency
f           = 0:150; % Used frequencies
sl_freq     = 12;    % Stimulus-locked frequency
sl_freq_i   = sl_freq + 1;    % We need to use the sl_freq +1 to get the correct index

tol         = 1.5;     % exclude frequencies within +/- tol of sl_freq
sl_drop     = f(mod(f, sl_freq) <= tol | mod(f, sl_freq) > sl_freq - tol);
   
% Exclude all frequencies that are close to a multiple of the
% line noise frequency (60 Hz)
ln_drop     = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lf_drop     = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);

keep_frequencies    = @(x) x(ab_i);
bb_frequencies      = f(ab_i);

% Define options for denoising
opt.resampling        = {'boot','boot'};
opt.pcselmethod       = 'snr';
opt.preprocessfun     = @(x) bbFilter(x, bb_frequencies);  % preprocess data with a high pass filter for broadband analysis
opt.npoolmethod       = {'snr','n',75};
opt.verbose           = true;

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,sl_freq_i); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,keep_frequencies,1000);  % function handle to compuite broadband with a sample rate of 1 kHz

% Get different epoch lengths and npcs to denoise with
epochDurs             = [1,3,6,12,24,36,72,inf];
npcs                  = [5,10:10:70];

% ------------------------------------------------------------------------
% ------------------ Load and denoise data per subject -------------------
% ------------------------------------------------------------------------

for whichSubject = subjects
    
    % ------------------ Load data and design ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    clear tmp;
    
    % ------------------ Make design matrix ------------------------------
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is left field
    design(conditions==7,3) = 1; % condition 7 is right field
    % condition 3 is blank
    
    
    %%
        % ------------------ Preprocess data ---------------------------------
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
    
    % ---- Define first epochs in order to remove later ------------------
    if removeFirstEpoch, badEpochs(1:6:end) = 1; end
    
    % -------------- Remove bad epochs and channels ----------------------
    sensorData      = sensorData(:,~badEpochs, ~badChannels);
    design          = design(~badEpochs,:);
    
    % ------------------ Permute sensorData for denoising ----------------
    sensorData = permute(sensorData, [3 1 2]);  
    
    
    %%
    
    
    
    % ------------------ Loop over Epoch length & nr PCs -----------------
    allResults = cell(length(epochDurs), length(npcs));
    for dur = 1:length(epochDurs) % iterate through epoch duration fo doing PCA 
        if isinf(epochDurs(dur))
            epochGroup = ones(size(sensorData,3),1);
        else
            epochGroup = megEpochGroup(~badEpochs,epochDurs(dur),0);
        end
        
        opt.epochgroup = epochGroup;
        clear results;
        fprintf('epochDur = %d; %d epochs\n', epochDurs(dur), length(epochGroup));
        
        for np = 1:length(npcs) % iterate through number of PCs projected out 
            opt.pcchoose   = -npcs(np);
            opt.npcs2try   = [];
            fprintf('\tnpcs = %d\n', npcs(np));
            
            if np == 1
                [results] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
                noisepooldef = results.noisepool;
            else
                [results] = denoisedata(design,sensorData,noisepooldef,evalfun,opt);
            end
            % convert function handle to string otherwise saved Matlab file
            % becomes orders of magnitude larger than necessary
            results.opt.preprocessfun = func2str(results.opt.preprocessfun);
            allResults{dur,np} = results; 
        end
    end
    
% ------------------------------------------------------------------------
% -------------------- Denoise and save the data -------------------------
% ------------------------------------------------------------------------
    
    fname = sprintf(fullfile(dfdRootPath,'analysis','data', 's0%d_denoisedData_varyEpochLength_NrPCs'),whichSubject);   
        
    parsave([fname '_bb.mat'], 'allResults', allResults, 'epochDurs', epochDurs, 'npcs', npcs);
    fprintf('Data saved:%s\n', fname);
    
end

  
end