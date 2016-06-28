%% s_ERF_analysis

% Script to plot the eye related function of the data

subjects = 6;

% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];%[0.05 30];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
use3Channels        = false;
removeFirstEpoch    = true;

%% Get frequencies to define stimulus locked and asynchronous broadband power
% Data are sampled at 1000 Hz and epoched in one-second bins. Hence
% frequencies are resolved in 1 Hz increments, 0:999 Hz
f           = 0:150;   % limit frequencies to [0 150] Hz
sl_freq     = 12;      % Stimulus-locked frequency
sl_freq_i   = sl_freq + 1;
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

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,sl_freq_i); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,keep_frequencies,1000);  % function handle to compuite broadband with a sample rate of 1 kHz

% Define options for denoising that are equal for each type of denoising
opt.resampling      = {'boot','boot'};
opt.pcselmethod     = 'snr';
opt.verbose         = true;

for whichSubject = subjects
    
    if whichSubject < 9,        dataChannels        = 1:157; % yokogawa MEG
    elseif whichSubject < 21,   dataChannels        = 1:204; % Elekta Neuromag
    else                        dataChannels        = 1:157; % Synthetic subject
    end

    % ------------------ Load data and design ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 'eye', 's%02d_ms.mat'),whichSubject)); ms = tmp.ms;
    
    conditions = [];
    for ii = 1:numel(ms)
        tmp = [ms{ii};ones(length(ms{ii})*ii,1)];
        conditions = vertcat([conditions;tmp]);
    end
    
    % ------------------ Make design matrix ------------------------------
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is left field
    design(conditions==7,3) = 1; % condition 7 is right field
    % condition 3 is blank
    
    % ------------- Combine channels if NeuroMag360 data -------------
    if whichSubject > 8
        optbb.npoolmethod = {'r2','n',100};
        optsl.npoolmethod = {'r2','n',100};
    end
    
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
    
end