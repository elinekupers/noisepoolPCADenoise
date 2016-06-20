function dfdDenoiseNPCvsNoisepool(whichSubjects)
% denoise as a function of denoising epoch duration and number of PCs
% removed
% 


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
f           = 0:150;   % limit frequencies to [0 150] Hz
sl_freq     = 12;      % Stimulus-locked frequency
sl_freq_i   = sl_freq + 1;
tol         = 1.5;     % exclude frequencies within +/- tol of sl_freq
sl_drop     = f(mod(f, sl_freq) <= tol | mod(f, sl_freq) > sl_freq - tol);
   
% Exclude all frequencies that are close to a multiple of the
% line noise frequency (60 Hz)
ln_drop     = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lf_drop = f(f<60);

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

% Denoise with exactly 10 PC regressors
opt.preprocessfun     =  @(x) bbFilter(x, bb_frequencies);
opt.verbose           = true;
npools                = [5,10:10:140];
npcs                  = [5,10:10:130];
allResults            = {};

for whichSubject = whichSubjects
    
    opt.preprocessfun     =  @(x) bbFilter(x, bb_frequencies);
    
    % ------------------ Load data and design ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    % ------------------ Make design matrix ------------------------------
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is left field
    design(conditions==7,3) = 1; % condition 7 is right field
    % condition 3 is blank
    
    % ------------------ Preprocess data ---------------------------------
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
    
        
    % ---- Define first epochs in order to remove later ------------------
    if removeFirstEpoch, badEpochs(1:6:end) = 1; end
    
    % ---- Remove bad channels and bad epochs from data and conditions ---
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design     = design(~badEpochs,:);
    
    % ------------------ Permute sensorData for denoising ----------------
    sensorData = permute(sensorData, [3 1 2]);  
    
    % Loop over the different eval functions
    for np = 1:length(npools)
         
        for nc = 1:length(npcs)
            fprintf('Denoising.. Using %d channels in noisepool and %d pcs to denoise \n',npools(np), npcs(nc))
            if npcs(nc)>npools(np), continue; end
            opt.npcs2try      = [];    
            opt.npoolmethod   = {'r2','n',npools(np)};
            opt.pcchoose      = -npcs(nc);
            [results] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
            
            % convert function handle to string otherwise saved Matlab file
            % becomes orders of magnitude larger than necessary
            results.opt.preprocessfun = func2str(results.opt.preprocessfun);
            
            allResults{np,nc} = results;          
        end
    end
    
    fname     = sprintf(fullfile(dfdRootPath,'analysis','data', ['s0%d_denoisedData_NCPSvsNoisePool']), whichSubject);
    parsave([fname '_bb.mat'], 'allResults', allResults)


end
