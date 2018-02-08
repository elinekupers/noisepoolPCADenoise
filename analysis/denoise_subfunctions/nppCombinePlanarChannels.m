function sensorDataCombined = nppCombinePlanarChannels(sensorData, dim)

if ~exist('dim', 'var'), 
    dim = find(size(sensorData) == 204);
end

% Check dimensions - we exect 204 channels to be combined to 102
assert(size(sensorData, dim) == 204);


if dim == 1, sensorData = sensorData'; end

sensorDataCombined = cat(3,sensorData(:,1:2:end), sensorData(:,2:2:end));
sensorDataCombined = nanmean(sensorDataCombined,3);

if dim == 1, sensorDataCombined = sensorDataCombined'; end

return
