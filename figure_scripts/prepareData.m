function [data,design,exampleIndex] = prepareData(dataDir,whichSubject)

% Load denoised data
load(sprintf(fullfile(dataDir, 's0%d_denoisedData_bb.mat'),whichSubject));
load(sprintf(fullfile(dataDir, 's0%d_conditions.mat'),whichSubject));
load(sprintf(fullfile(dataDir, 's0%d_sensorData.mat'),whichSubject));
                                         
%% Spectrum before and after denoising

% Because some channels were removed prior to the denoising (due to
% excessive noise), there are fewer channels in the denoised data than in
% the original data. We need to identify the correspondence.
channelNumbers = find(~badChannels);
exampleChannel = 42;
exampleIndex   = find(channelNumbers == exampleChannel);

sensorData = permute(sensorData, [3 1 2]);
% time domain data before and after denoising
data = {sensorData,denoisedts{1}}; %#ok<USENS>

% Compute design matrix
design = zeros(size(conditions,1),3);
design(conditions == 1,1) = 1; % Full
design(conditions == 5,2) = 1; % Right
design(conditions == 7,3) = 1; % Left

design = {design,design(~badEpochs,:)};
