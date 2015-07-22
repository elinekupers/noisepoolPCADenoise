function dfdDenoiseDifferentNPCsNoisePools(whichSubjects)
% denoise as a function of denoising epoch duration and number of PCs
% removed
% you can modify the script take parameter k for session number



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
[ab_f, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);

keep_frequencies    = @(x) x(ab_f);

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,sl_freq); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,keep_frequencies,1000);  % function handle to compuite broadband with a sample rate of 1 kHz

% Define options for denoising that are equal for each type of denoising
opt.resampling      = {'boot','boot'};
opt.pcselmethod     = 'snr';

% Denoise with exactly 10 PC regressors
opt.pcchoose          = -10;   % Take 10 PCs
opt.npcs2try          = [];    
optbb                 = opt;
optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
npools                = [5,10:10:140];
npcs                  = [5,10:10:130];
allResults            = [];

for whichSubject = whichSubjects
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
    
    % Loop over the different eval functions
%     for kk = 1:length(evalfun), evalfunstr{kk} = func2str(evalfun{kk}); end %#ok<AGROW>

    for np = 1:length(npools)
        for nc = 1:length(npcs)
            if npcs(nc)>npools(np), continue; end
            
            opt.npcs2try      = [];    
            opt.npoolmethod   = {'snr','n',npools(np)};
            opt.pcchoose      = -npcs(nc);
            [results,evalout] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
            allResults(np,nc) = results; %#ok<AGROW>
        end
    end
    
    postFix   = evalfunstr{kk};
    fname     = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data', ['s0%d_denoisedData_NCPSvsNoisePool_' postFix]), whichSubject);
    parsave([fname '_bb.mat'], 'results', results, 'evalout', evalout, 'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)


end
