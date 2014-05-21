function okEpochs = megIdenitfyBadEpochs(sensorDataIn, threshold)
% sensorDataOut = eegRemoveBadEpochs(sensorDataIn, threshold)
%
% Remove epochs with when the number of bad channels in that epoch exceeds threshold
if notDefined('threshold'), threshold = 0.5; end

% sensor data is channels x time x epoch
badData = squeeze(isnan(sensorDataIn(:,1,:)));
% average across channels
okEpochs = mean(badData)' < threshold;
