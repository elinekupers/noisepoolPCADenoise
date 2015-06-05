function sensorDataOut = dfdChannelRepair(sensorDataIn, outliers, method, sensorPositions)
% Replace bad data with an interpolation from good channels
%  
%sensorDataOut = dfdChannelRepair(sensorDataIn, outliers, method)
%
% INPUTS
%   sensorDataIn: data array, time points x epochs x channels
%   outliers:     binary matrix,epochs x channels. 1 = bad, 0 = good.
%   method:       interpolation method. Currently must be 'nearest'. In
%                   principle, could allow for splines, average, laplacian,
%                   etc.
% OUTPUTS
%   sensorDataOut: data array, same size as sensorDataIn

if notDefined('sensorPositions'), 
    hdr = load('meg160_example_hdr.mat'); hdr = hdr.hdr;
    net.xyz = hdr.grad.chanpos; 
end

% get data sizes
nEpochs   = size(sensorDataIn,2);
nChannels = size(sensorDataIn,3);

% compute distance matrix, which will be used for weighting channels for
% interpolation
connectivityMatrix = eye(nChannels);
distances          = squareform(pdist(net.xyz), 'tomatrix');  
epochsAllBad       = all(outliers, 2);
epochsAllGood      = all(~outliers, 2); 
epochsToInterpolate = find(~epochsAllBad & ~epochsAllGood); % these epochs have no bad data, so do not interpolate

% Initialize the sensorDataOut matrix
sensorDataOut         = zeros(size(sensorDataIn));

 % These epochs have no usable data.
sensorDataOut(:,epochsAllBad,:) = NaN;

% These epochs have no bad data, so do not interpolate
sensorDataOut(:,epochsAllGood,:) = sensorDataIn(:,epochsAllGood,:);

% Interpolate the rest
for ii = 1:length(epochsToInterpolate)
    thisEpoch = epochsToInterpolate(ii);
    badChannels = outliers(thisEpoch,:);
    goodChannels = ~badChannels;
    
    weightMatrix = connectivityMatrix;
    weightMatrix(:,badChannels) = 0;
    weightMatrix(badChannels, goodChannels) = 1./distances(badChannels, goodChannels);
    weightMatrix = bsxfun(@rdivide, weightMatrix, sum(weightMatrix,2));
    
    thisdata = weightMatrix*permute(sensorDataIn(:,thisEpoch,:), [3 1 2]);
    sensorDataOut(:,thisEpoch,:) = thisdata';
end
