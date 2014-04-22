function badChannels = eegIdenitfyBadEpochs(sensorData, threshold)
% sensor data is channels x time x epoch
%
%  badChannels = eegIdenitfyBadChannels(sensorData)

if ~iscell(sensorData), sensorData = {sensorData}; end


badChannels = zeros(nmbChannels, 1);

for ii = 1:length(sensorData)
    nmbEpochs = size(sensorData{ii}, 1);
    
    % if an epoch is OK, then all time points should be finite. if the
    % epoch is not OK, then all time points should be NaN. so just check
    % the first time point for each epoch of each channel
    tmp  = squeeze(sensorData{ii}(:,1,:));
    okepochs    = isfinite(tmp);
    badChannels = badChannels | mean(okepochs,2) < threshold;
end
return
