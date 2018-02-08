function sensorDataOut = nppChannelRepair(sensorDataIn, outliers, method, sensorPositions, verbose)
% Replace bad data with an interpolation from good channels
%  
%sensorDataOut = nppChannelRepair(sensorDataIn, outliers, method)
%
% INPUTS
%   sensorDataIn: data array, time points x epochs x channels
%   outliers:     binary matrix,epochs x channels. 1 = bad, 0 = good.
%   method:       interpolation method. Currently must be 'nearest'. In
%                   principle, could allow for splines, average, laplacian,
%                   etc.
% OUTPUTS
%   sensorDataOut: data array, same size as sensorDataIn
if notDefined('verbose'); verbose = false; end
if notDefined('sensorPositions') 
    hdr = load('meg160_example_hdr.mat'); hdr = hdr.hdr;
    net.xyz = hdr.grad.chanpos;
elseif strcmp(sensorPositions, 'meg160_example_hdr');
    hdr = load('meg160_example_hdr.mat'); hdr = hdr.hdr;
    net.xyz = hdr.grad.chanpos; 
elseif strcmp(sensorPositions, 'neuromag360xyz')
    load('neuromag360_sample_cfg_combined');    
    channels = cfg.cfg.channel;
    hdr = load('neuromag360_sample_hdr.mat'); hdr = hdr.hdr;
    idx = NaN(1,length(channels));
    for ii = 1:length(channels)
       idx(ii) = find(strcmp(channels{ii}, hdr.grad.label));
    end
    
    net.xyz = hdr.grad.chanpos(idx,:);
elseif strcmp(sensorPositions, 'yokogawa_con_example_hdr')
    hdr = load('yokogawa_con_example_hdr.mat'); hdr = hdr.hdr;
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

% for planar gradiometers, we may have paired gradiometers with 0 distance
% between them. this will cause a problem when normalizing weights by
% distance. so replace 0 with a small number
mn = min(distances(distances(:)>0));
distances(distances==0) = mn/2;
distances = distances .* (1-eye(nChannels));

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
