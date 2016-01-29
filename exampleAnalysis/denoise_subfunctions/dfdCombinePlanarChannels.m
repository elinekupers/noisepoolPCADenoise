function [results, evalout] = dfdCombinePlanarChannels(whichSubject, denoisedts, design, badEpochs, sl_freq_i, keep_frequencies)

% Get fieldtrip cfg structure from raw data
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data'); 

% Load data
dataRaw = load(sprintf(fullfile(dataDir, 's%02d_raw.mat'),whichSubject));
conditions = load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));

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
dataRaw.trial = epoch2trial(dataRaw.trial, denoisedts{1}(:,:,epochs));

% Combine planar gradiometers
dataRaw.cfg.demean  = 'no';
dataRaw.cfg.trials  = 'all';
dataCombined    = ft_combineplanar(dataRaw.cfg, dataRaw);

% Transform back to one array: channels x time x epochs
denoisedts = catcell(3, dataCombined.trial);



% Run algorithm again without denoising to get results back in the right
% format
opt.preprocessfun = [];
opt.npcs2try    = 0;
opt.pcchoose    = 0;
evokedfun       = @(x)getstimlocked(x,sl_freq_i);
evalfun         = @(x)getbroadband(x,keep_frequencies,1000);

[results,evalout] = denoisedata(design,denoisedts,evokedfun,evalfun,opt);

return
