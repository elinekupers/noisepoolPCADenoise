function [sensorDataOut,okEpochsOut] = eegRemoveBadEpochs(sensorDataIn, threshold)
% sensorDataOut = eegRemoveBadEpochs(sensorDataIn, threshold)
%
% Remove epochs with when the number of bad channels in that epoch exceeds threshold
if notDefined('threshold'), threshold = 0.5; end

sensorDataOut = cell(size(sensorDataIn));
okEpochsOut   = cell(size(sensorDataIn));
for ii = 1:numel(sensorDataIn)
    
    % sensor data is channels x time x epoch
    sz = size(sensorDataIn{ii});
    badData = squeeze(isnan(sensorDataIn{ii}(:,1,:)));
    okEpochs = mean(badData) < threshold;
    sensorDataOut{ii} = sensorDataIn{ii}(:,:,okEpochs);
    okEpochsOut{ii}   = okEpochs;
    fprintf('cond %d: okEpochs: %d / %d\n', ii, sum(okEpochs), length(okEpochs));
end

return