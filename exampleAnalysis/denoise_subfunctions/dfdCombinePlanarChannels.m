function sensorDataCombined = dfdCombinePlanarChannels(whichSubject, sensorData, badEpochs)

% Get fieldtrip cfg structure from raw data
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data'); 

% Load data
dataRaw = load(sprintf(fullfile(dataDir, 's%02d_raw.mat'),whichSubject));

% Combine different fields in struct (TO DO: make general statement)
combined = [dataRaw.data_A_all, dataRaw.data_R_all, dataRaw.data_L_all];

dataRaw = combined(1);

% Get epochs, and remove the ones that were defined as bad
epochs = 1:(12*size(combined(1).trial,2) + 12*size(combined(2).trial,2) + 12*size(combined(3).trial,2));
epochs = epochs(~badEpochs);

% Use knowledge of number of epochs to shorten number of cells in
% dataRaw.trial
dataRaw.trial = repmat({zeros(204,1000)},size(epochs,2),1);

% Put the denoised data (epoched into 1000 ms epochs) into fieldtrip struct
dataRaw.trial = epoch2trial(dataRaw.trial, sensorData);

% Combine planar gradiometers
dataRaw.cfg.demean  = 'no';
dataRaw.cfg.trials  = 'all';
dataCombined    = ft_combineplanar(dataRaw.cfg, dataRaw);

% Transform back to one array: channels x time x epochs
sensorDataCombined = catcell(3, dataCombined.trial);

return
