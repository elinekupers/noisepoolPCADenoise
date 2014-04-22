function badChannels = eegIdenitfyBadChannels(sensorData, thresh)
% sensor data is channels x time x epoch
%
%  badChannels = eegIdenitfyBadChannels(sensorData)

if ~iscell(sensorData), sensorData = {sensorData}; end

nmbChannels = size(sensorData{1}, 1);
badChannels = zeros(nmbChannels, 1);

for ii = 1:length(sensorData)
    
    % if an epoch is OK, then all time points should be finite. if the
    % epoch is not OK, then all time points should be NaN. so just check
    % the first time point for each epoch of each channel
    tmp  = squeeze(sensorData{ii}(:,1,:));
    okepochs    = isfinite(tmp);
    badChannels = badChannels | mean(okepochs,2) < thresh;
end
return
