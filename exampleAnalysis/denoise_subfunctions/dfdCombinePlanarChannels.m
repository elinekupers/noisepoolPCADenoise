function sensorDataCombined = dfdCombinePlanarChannels(whichSubject, sensorData)

% Permute sensorData
sensorData = permute(sensorData, [3 1 2]);

% Get fieldtrip cfg structure from raw data
dataDir         = fullfile(dfdRootPath, 'exampleAnalysis', 'data'); 

% Load data
dataRaw = load(sprintf(fullfile(dataDir, 's%02d_raw.mat'),whichSubject));

% Combine different fields in struct (TO DO: make general statement)
% if size(dataRaw.data_A_all.trial) ~= size(dataRaw.data_R_all.trial)
%     

try 
    combined = cat(1,[dataRaw.data_A_all, dataRaw.data_R_all, dataRaw.data_L_all]);
catch 
    warning('(dfdCombinePlanarChannels): Sample info was not present in one of the data set and will be added');
    
    dataRaw.data_R_all.sampleinfo = dataRaw.data_A_all.sampleinfo;
    dataRaw.data_L_all.sampleinfo = dataRaw.data_A_all.sampleinfo;
    
    combined = cat(1,[dataRaw.data_A_all, dataRaw.data_R_all, dataRaw.data_L_all]);

end

dataRaw = combined(1);

% Get epochs, and remove the ones that were defined as bad
epochs = 1:(12*size(combined(1).trial,2) + 12*size(combined(2).trial,2) + 12*size(combined(3).trial,2));
% Use knowledge of number of epochs to shorten number of cells in
% dataRaw.trial
dataRaw.trial = repmat({zeros(204,1000)},size(epochs,2),1);
dataRaw.time  = repmat({1:1:1000}, size(epochs,2),1);


% Put the denoised data (epoched into 1000 ms epochs) into fieldtrip struct
dataRaw.trial = epoch2trial(dataRaw.trial, sensorData);

% Combine planar gradiometers
dataRaw.cfg.demean  = 'no';
dataRaw.cfg.trials  = 'all';
dataCombined    = ft_combineplanar(dataRaw.cfg.trials, dataRaw);

% Transform back to one array: channels x time x epochs
sensorDataCombined = catcell(3, dataCombined.trial);

sensorDataCombined = permute(sensorDataCombined, [2 3 1]);

return
