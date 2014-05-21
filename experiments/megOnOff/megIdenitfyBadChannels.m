function badChannels = megIdenitfyBadChannels(sensorData, thresh)
% sensor data is channels x time x epoch
%
%  badChannels = eegIdenitfyBadChannels(sensorData)

% if an epoch is OK, then all time points should be finite. if the
% epoch is not OK, then all time points should be NaN. so just check
% the first time point for each epoch of each channel
tmp  = squeeze(sensorData(:,1,:));
okepochs    = isfinite(tmp);
badChannels = mean(okepochs,2) < thresh;
