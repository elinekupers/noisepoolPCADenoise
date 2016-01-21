function [results, evalout] = dfdCombinePlanarChannels(whichSubject, denoisedts, badEpochs)

% Get fieldtrip cfg structure from raw data
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data'); 

% Load data
dataRaw = load(sprintf(fullfile(dataDir, 's%02d_raw.mat'),whichSubject));
conditions = load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));

% Combine different fields in struct (TO DO: make general statement)
combined = [dataRaw.data_A_all, dataRaw.data_R_all, dataRaw.data_L_all];

% Get epochs, and remove the ones that were defined as bad
epochs = 1:(12*size(combined(1).trial,2) + 12*size(combined(2).trial,2) + 12*size(combined(3).trial,2));
epochs = epochs(~badEpochs);

epochs = 

% Put the denoised data (epoched into 1000 ms epochs) into fieldtrip struct
dataRaw.cfg.trial = epoch2trial(epochs, denoisedts{1});

% Combine planar gradiometers
dataRaw.cfg.demean  = 'no';
dataCombined    = ft_combineplanar(dataRaw.cfg, dataRaw);

% Transform back to epochs
sensorData = trial2epoch(dataRaw.trial);

% Run algorithm again without denoising to get results back in the right
% format
opt.npcs2try    = 0;
opt.pcchoose    = 0;
evokedfun       = @(x)getstimlocked(x,opt.freq);
evalfun         = @(x)getbroadband(x,opt.freq);

[results,evalout] = denoisedata(conditions,sensorData,evokedfun,evalfun,opt);

return
