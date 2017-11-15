function dfdMakeSyntheticDataSet_pinknoise()

% Function to create a synthetic data set 
%
% dfdMakeSyntheticDataSet()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)
% 
% The synthetic data set was comprised of a mixture of 4 components in each
% of 157 sensors across 30 epochs (1 s each, with ms sampling). The 4
% components were stimulus locked, broadband, background uncorrelated
% noise, and background correlated noise. We will save the synthetic data
% set (sensorData and conditions mat file) as subject 99.

clear; close all;

%% --------------------- Specify simulation parameters --------------------

%   Epochs will be one second long, sampled at 1000 Hz. The stimuli are
%   chosen to match those in the experiment: Both-field stimulation
%   (Both=1), left hemifield stimulation (L=5), right hemifield stimulation
%   (R=7), and blank (Blank=3). The numbers 1, 5, 7, and 3 are essentially
%   arbitrary. (They happen to be the trigger codes used during the actual
%   experiment for full field, left, right, and blank conditions). For
%   conistency, we use these same values here. The sequence is similar to
%   the experiment but truncated, conistent of 6-second blocks of
%   stimuli alternating with 6-second blanks. Because the data are epoched
%   into one-second periods for analysis, we list the conditions as
%   FFFFFF-BBBBBB-LLLLLL-BBBBBB-RRRRRR-BBBBBB

fs                            = 1000; % Sampling frequency (Hz)
epochLength                   = 1;    % Seconds
blockLength                   = 6;    % Seconds
numChannels                   = 157;
numConditions                 = 4; 	  % Both-field, left, right, blank
conditionCodes                = [1 5 7 3]; % Arbitrary codes assigned to conditions
conditionNames                = {'both-field' 'left-field' 'right-field' 'blank'};
blockOrder                    = [1 4 2 4 3 4]; % Both/blank, left/blank, right/blank
numRepeats                    = 1;
stimulusLockedFrequency       = 12;
numNoiseBasis                 = 10;    % Number of independent bases for correlated noise
responseAmp.uncorrelatedNoise = 1;
responseAmp.stimulusLocked    = 2;
responseAmp.broadband         = 2;
responseAmp.correlatedNoise   = 5;

% Where to save data?
dataDir                       = fullfile(dfdRootPath, 'analysis', 'data');    

% ---------------------------- Derived ------------------------------------
numBlocks           = length(blockOrder) * numRepeats;
epochsPerBlock      = blockLength / epochLength;
numEpochs           = numBlocks * epochsPerBlock;
blockOrderAll       = repmat(blockOrder, 1, numRepeats);
samplesPerEpoch     = epochLength * fs;
epochsOrder         = reshape(ones(epochsPerBlock,1)*blockOrderAll, [], 1);
conditions          = conditionCodes(epochsOrder)';


%% Group channels to noisepool, visually responding left, visually responding right 

% Load Yokogawa MEG header for channel positions
load('meg160_example_hdr');
xyz = hdr.grad.chanpos;

% Find the 75 channels in front (highest x value): these will be noise pool
[~, sortY] = sort(xyz(:,1), 'descend');
sensors.noisepool = sortY(1:75); 

% Among the 82 non-noise pool, find the 41 channels that are most left and
% the 41 that are most right
nonNoisepool = sortY(76:end);
[~, sortX] = sort(xyz(nonNoisepool,2), 'descend');
sensors.right  = nonNoisepool(sortX(1:41));
sensors.left = nonNoisepool(sortX(42:end));

% Create map
map = zeros(1, 157);
map(sensors.noisepool) = 1;
map(sensors.left) = 2;
map(sensors.right) = 3;

% Plot it
megPlotMap(map)

%% -------------------- Generate noiseless time series --------------------

% Predefine sensorData array
sensorData = zeros(samplesPerEpoch, numEpochs, numChannels);

% Define time points
t = (1:1000)/fs;

% ------------------- COMPONENT 1: STIMULUS LOCKED ------------------------
% Make an example stimulus-locked response
stimulusLocked =  responseAmp.stimulusLocked * 2 * ...
    sin(2*pi*t*stimulusLockedFrequency);

% For both field epochs, fill the left and right sensor data with stimulus
% locked time series
sz = size(sensorData(:,conditions==1, [sensors.left sensors.right]));
sensorData(:,conditions==1, [sensors.left sensors.right]) = ...
    repmat(stimulusLocked(:), [1 sz(2) sz(3)]);

% For left field epochs, fill the left sensor data with stimulus locked time series
sz = size(sensorData(:,conditions==5, [sensors.left]));
sensorData(:,conditions==5, [sensors.left]) = ...
    repmat(stimulusLocked(:), [1 sz(2) sz(3)]);

% For right field epochs, fill the left sensor data with stimulus locked time series
sz = size(sensorData(:,conditions==7, [sensors.right]));
sensorData(:,conditions==7, [sensors.right]) = ...
    repmat(stimulusLocked(:), [1 sz(2) sz(3)]);

% ------------------ COMPONENT 2: BROADBAND -------------------------------
% For both field epochs, fill the left and right sensor data with broadband time series
sz = size(sensorData(:,conditions==1, [sensors.left sensors.right]));
sensorData(:,conditions==1, [sensors.left sensors.right]) = ...
    sensorData(:,conditions==1, [sensors.left sensors.right]) + ...
    responseAmp.broadband * zscore(generatepinknoisedata(sz));

% For left field epochs, fill the left sensor data with broadband time series
sz = size(sensorData(:,conditions==1, [sensors.left]));
sensorData(:,conditions==5, [sensors.left]) = ...
    sensorData(:,conditions==5, [sensors.left]) + ...
    responseAmp.broadband * zscore(generatepinknoisedata(sz));

% For right field epochs, fill the right sensor data with broadband time series
sz = size(sensorData(:,conditions==7, [sensors.right]));
sensorData(:,conditions==7, [sensors.right]) = ...
    sensorData(:,conditions==7, [sensors.right]) + ...
    responseAmp.broadband * zscore(generatepinknoisedata(sz));

% ------------------- COMPONENT 3: UNCORRELATED NOISE  --------------------
% Add Gaussian white noise to all sensors, all conditions
sz = [samplesPerEpoch, numEpochs, numChannels];
sensorData = sensorData + ...
    responseAmp.uncorrelatedNoise * zscore(generatepinknoisedata(sz));


% ------------------- COMPONENT 4: CORRELATED NOISE  ----------------------
% Create basis functions that are correlated
%correlatedBasis = randn(samplesPerEpoch * numEpochs, numNoiseBasis);
sz = [samplesPerEpoch, numEpochs, numNoiseBasis];
correlatedBasis = zscore(generatepinknoisedata(sz));
correlatedBasis = reshape(correlatedBasis, [], numNoiseBasis);
mixingMatrix    = randn(numNoiseBasis, numChannels);
mixingMatrix    = bsxfun(@rdivide, mixingMatrix, sqrt(sum(mixingMatrix.^2)));
correlatedNoise = correlatedBasis * mixingMatrix;
correlatedNoise = reshape(correlatedNoise, samplesPerEpoch, numEpochs, numChannels);

sensorData       = sensorData + correlatedNoise * responseAmp.correlatedNoise;
    
%% ------------------------- Save it --------------------------------------
save(fullfile(dataDir, 's99_sensorData'), 'sensorData');
save(fullfile(dataDir, 's99_conditions'), 'conditions');



end

function pinkdata = generatepinknoisedata(sz)

    alpha = 1; % = 1/f (pink noise)
    for chan = 1:sz(3)
        pinkdata(:,:,chan) = generatepinknoise1D(sz(1),alpha,sz(2),1);
    end

end


