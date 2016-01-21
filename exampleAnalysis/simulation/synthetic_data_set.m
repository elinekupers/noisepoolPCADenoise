% Create synthetic data set for validation of MEG denoisedata pipeline

%% specify simulation parameters

fs = 1000;           % sampling frequency, Hz
epoch_length    = 1; % seconds
block_length    = 6; % seconds
num_channels    = 157;
num_conditions  = 4; % full-field, left, right, blank
condition_codes = [1 5 7 3]; % arbitrary codes assigned to conditions
block_order     = [1 4 2 4 3 4]; % full/blank, left/blank, right/blank
num_repeats     = 1;
stimulus_locked_frequency = 12;
num_noise_basis                 = 10; % number of independent bases for correlated noise
response_amp.uncorrelated_noise = 1;
response_amp.stimulus_locked    = 2;
response_amp.broadband          = 2;
response_amp.correlated_noise   = 5;

% --------- derived -------------------------
num_blocks        = length(block_order) * num_repeats;
epochs_per_block  = block_length / epoch_length;
num_epochs        = num_blocks * epochs_per_block;
condition_codes   = [1 5 7 3];
condition_names   = {'full-field' 'left-field' 'right-field' 'blank'};
block_order_all   = repmat(block_order, 1, num_repeats);
samples_per_epoch = epoch_length * fs;
epochs_order      = reshape(ones(epochs_per_block,1)*block_order_all, [], 1);

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

conditions = condition_codes(epochs_order)';
save('s99_conditions', 'conditions');

%% Group channels to noisepool, respond left, respond right
load('meg160xyz');

% find the 75 channels in front (highest x value): these will be noise pool
[~, sort_y] = sort(xyz(:,1), 'descend');
channels.noisepool = sort_y(1:75); 

% among the 82 non-noise pool, find the 41 channels that are most left and
% the 41 that are most right
nonnoisepool = sort_y(76:end);
[~, sort_x] = sort(xyz(nonnoisepool,2), 'descend');
channels.left  = nonnoisepool(sort_x(1:41));
channels.right = nonnoisepool(sort_x(42:end));

% % check it
% map = zeros(1, 157);
% map(noisepool) = 1;
% map(left) = 2;
% map(right) = 3;
% 
% ft_plotOnMesh(map, [],[],[],'interpolation', 'nearest')

%% generate noiseless time series
sensorData = zeros(samples_per_epoch, num_epochs, num_channels);

t = (1:1000)/fs;

% ------------------------STIMULUS LOCKED --------------------------------
%
% make an example stimulus-locked response
stimulus_locked =  response_amp.stimulus_locked * 2 * ...
    sin(2*pi*t*stimulus_locked_frequency);

% for full field epochs, fill the left and right sensor data with stimulus
% locked time series
sz = size(sensorData(:,conditions==1, [channels.left channels.right]));
sensorData(:,conditions==1, [channels.left channels.right]) = ...
    repmat(stimulus_locked(:), [1 sz(2) sz(3)]);

% for left field epochs, fill the left sensor data with stimulus locked time series
sz = size(sensorData(:,conditions==5, [channels.left]));
sensorData(:,conditions==5, [channels.left]) = ...
    repmat(stimulus_locked(:), [1 sz(2) sz(3)]);

% for right field epochs, fill the left sensor data with stimulus locked time series
sz = size(sensorData(:,conditions==7, [channels.right]));
sensorData(:,conditions==7, [channels.right]) = ...
    repmat(stimulus_locked(:), [1 sz(2) sz(3)]);

% ------------------------BROADBAND  --------------------------------
% 
% for full field epochs, fill the left and right sensor data with broadband time series
sz = size(sensorData(:,conditions==1, [channels.left channels.right]));
sensorData(:,conditions==1, [channels.left channels.right]) = ...
    sensorData(:,conditions==1, [channels.left channels.right]) + ...
    response_amp.broadband * zscore(randn(sz));

% for left field epochs, fill the left sensor data with broadband time series
sz = size(sensorData(:,conditions==1, [channels.left]));
sensorData(:,conditions==5, [channels.left]) = ...
    sensorData(:,conditions==5, [channels.left]) + ...
    response_amp.broadband * zscore(randn(sz));

% for right field epochs, fill the right sensor data with broadband time series
sz = size(sensorData(:,conditions==7, [channels.right]));
sensorData(:,conditions==7, [channels.right]) = ...
    sensorData(:,conditions==7, [channels.right]) + ...
    response_amp.broadband * zscore(randn(sz));

% ------------------------Uncorrelated Noise  ----------------------------
% add white noise to all channels, all conditions
sensorData = sensorData + ...
    response_amp.uncorrelated_noise * zscore(randn(samples_per_epoch, num_epochs, num_channels));


% ------------------------Correlated Noise  ----------------------------
% 

correlated_basis = randn(samples_per_epoch * num_epochs, num_noise_basis);
mixing_matrix    = randn(num_noise_basis, num_channels);
mixing_matrix    = bsxfun(@rdivide, mixing_matrix, sqrt(sum(mixing_matrix.^2)));
correlated_noise = correlated_basis * mixing_matrix;
correlated_noise = reshape(correlated_noise, samples_per_epoch, num_epochs, num_channels);

sensorData       = sensorData + correlated_noise * response_amp.correlated_noise;
    
% ------------------------ Save it -------------------------------------
save('s99_sensorData', 'sensorData');

