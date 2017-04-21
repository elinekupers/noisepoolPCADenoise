% Create synthetic data set for validation of MEG denoisedata pipeline

%% specify simulation parameters

fs = 1000;           % sampling frequency, Hz
epochLength    = 1; % seconds
blockLength    = 6; % seconds
numChannels    = 157;
numConditions  = 4; % full-field, left, right, blank
conditionCodes = [1 5 7 3]; % arbitrary codes assigned to conditions
conditionNames   = {'full-field' 'left-field' 'right-field' 'blank'};
blockOrder     = [1 4 2 4 3 4]; % full/blank, left/blank, right/blank
numRepeats     = 1;
stimulusLockedFrequency = 12;
numNoiseBasis                 = 10; % number of independent bases for correlated noise
responseAmp.uncorrelatedNoise = 1;
responseAmp.stimulusLocked    = 2;
responseAmp.broadband         = 2;
responseAmp.correlatedNoise   = 5;

% --------- derived -------------------------
numBlocks        = length(blockOrder) * numRepeats;
epochsPerBlock  = blockLength / epochLength;
numEpochs        = numBlocks * epochsPerBlock;
blockOrderAll   = repmat(blockOrder, 1, numRepeats);
samplesPerEpoch = epochLength * fs;
epochsOrder      = reshape(ones(epochsPerBlock,1)*blockOrderAll, [], 1);

%   Epochs will be one second long, sampled at 1000 Hz. The stimuli are
%   chosen to match those in the example experiment: full-field stimulation
%   (F=1), left hemifield stimulation (L=5), right hemifield stimulation
%   (R=7), and blank (B=3). THe numbers 1, 5, 7, and 3 are essentially
%   arbitrary. (They happen to be the trigger codes used during the actual
%   experiment for full field, left, right, and blank conditions). For
%   conistency, we use these same values here. THe sequence is similar to
%   the example experiment but truncated, conistent of 6-second blocks of
%   stimuli alternating with 6-second blanks. Because the data are epoched
%   into one-second periods for analysis, we list the conditions as
%   FFFFFF-BBBBBB-LLLLLL-BBBBBB-RRRRRR-BBBBBB

conditions = conditionCodes(epochsOrder)';


%% Group channels to noisepool, respond left, respond right
load('meg160_example_hdr');
xyz = hdr.grad.chanpos;

% find the 75 channels in front (highest x value): these will be noise pool
[~, sortY] = sort(xyz(:,1), 'descend');
channels.noisepool = sortY(1:75); 

% among the 82 non-noise pool, find the 41 channels that are most left and
% the 41 that are most right
nonnoisepool = sortY(76:end);
[~, sortX] = sort(xyz(nonnoisepool,2), 'descend');
channels.right  = nonnoisepool(sortX(1:41));
channels.left = nonnoisepool(sortX(42:end));

% check it
map = zeros(1, 157);
map(channels.noisepool) = 1;
map(channels.left) = 2;
map(channels.right) = 3;

megPlotMap(map)

%% generate noiseless time series
sensorData = zeros(samplesPerEpoch, numEpochs, numChannels);

t = (1:1000)/fs;

% ------------------------STIMULUS LOCKED --------------------------------
%
% make an example stimulus-locked response
stimulusLocked =  responseAmp.stimulusLocked * 2 * ...
    sin(2*pi*t*stimulusLockedFrequency);

% for full field epochs, fill the left and right sensor data with stimulus
% locked time series
sz = size(sensorData(:,conditions==1, [channels.left channels.right]));
sensorData(:,conditions==1, [channels.left channels.right]) = ...
    repmat(stimulusLocked(:), [1 sz(2) sz(3)]);

% for left field epochs, fill the left sensor data with stimulus locked time series
sz = size(sensorData(:,conditions==5, [channels.left]));
sensorData(:,conditions==5, [channels.left]) = ...
    repmat(stimulusLocked(:), [1 sz(2) sz(3)]);

% for right field epochs, fill the left sensor data with stimulus locked time series
sz = size(sensorData(:,conditions==7, [channels.right]));
sensorData(:,conditions==7, [channels.right]) = ...
    repmat(stimulusLocked(:), [1 sz(2) sz(3)]);

% ------------------------BROADBAND  --------------------------------
% 
% for full field epochs, fill the left and right sensor data with broadband time series
sz = size(sensorData(:,conditions==1, [channels.left channels.right]));
sensorData(:,conditions==1, [channels.left channels.right]) = ...
    sensorData(:,conditions==1, [channels.left channels.right]) + ...
    responseAmp.broadband * zscore(randn(sz));

% for left field epochs, fill the left sensor data with broadband time series
sz = size(sensorData(:,conditions==1, [channels.left]));
sensorData(:,conditions==5, [channels.left]) = ...
    sensorData(:,conditions==5, [channels.left]) + ...
    responseAmp.broadband * zscore(randn(sz));

% for right field epochs, fill the right sensor data with broadband time series
sz = size(sensorData(:,conditions==7, [channels.right]));
sensorData(:,conditions==7, [channels.right]) = ...
    sensorData(:,conditions==7, [channels.right]) + ...
    responseAmp.broadband * zscore(randn(sz));

% ------------------------Uncorrelated Noise  ----------------------------
% add white noise to all channels, all conditions
sensorData = sensorData + ...
    responseAmp.uncorrelatedNoise * zscore(randn(samplesPerEpoch, numEpochs, numChannels));


% ------------------------Correlated Noise  ----------------------------
% 

correlatedBasis = randn(samplesPerEpoch * numEpochs, numNoiseBasis);
mixingMatrix    = randn(numNoiseBasis, numChannels);
mixingMatrix    = bsxfun(@rdivide, mixingMatrix, sqrt(sum(mixingMatrix.^2)));
correlatedNoise = correlatedBasis * mixingMatrix;
correlatedNoise = reshape(correlatedNoise, samplesPerEpoch, numEpochs, numChannels);

sensorData       = sensorData + correlatedNoise * responseAmp.correlatedNoise;
    
% ------------------------ Save it -------------------------------------
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
save(fullfile(dataDir, 's99_sensorData'), 'sensorData');
save(fullfile(dataDir, 's99_conditions'), 'conditions');
