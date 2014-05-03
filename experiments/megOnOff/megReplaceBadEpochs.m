function sensorDataOut = megReplaceBadEpochs(sensorDataIn,net,montageType,dist)
% sensorDataOut = eegReplaceBadChannels(sensorDataIn, montageType)


% Replace channels with mean of neighboring channels
if notDefined('montageType'), montageType = 'egi'; end
if notDefined('dist'), dist = 3; end
%if ~exist('defaultSphereNet.mat', 'file'), addpath('mrCurrent'); end

% sensor data is channels x time x epoch
sz = size(sensorDataIn);

switch lower(montageType)
    case {'128', 'egi', 'hd'}
        % we expect 157 channels
        assert(sz(1) == 157);
        
        % load the xyz positions of the sensor
        %net=load('defaultSphereNet.mat');
        
        % calculate a distance matrix
        distances = squareform(pdist(net.xyz), 'tomatrix');
        
        % for each channel, find neighbors within X mm
        nearest = distances < dist & distances > 0;
        
        % the sum of the columns should be one
        nearest = bsxfun(@rdivide, nearest, sum(nearest));
        
        % eliminate any NaNs (channels with no neighbors)
        nearest(isnan(nearest)) = 0;
        
        % reshape 3d sensor data into matrix (2d) in order to multiply
        sensorDataIn = reshape(sensorDataIn,  sz(1), []);
        
        % Initialize out data
        sensorDataOut = sensorDataIn;
        
        
        % Replace bad channels with the mean of the neighbors
        
        for jj = 1:157
            inds = find(nearest(:, jj));
            %if jj <= 3, disp(length(inds)); end
            badPoints = isnan(sensorDataOut(jj,:));
            tmp = sensorDataIn(inds,badPoints);
            tmp(isnan(tmp)) = 0;
            sensorDataOut(jj,badPoints) = (tmp' * nearest(inds, jj))';
        end
        % reshape 2d sensorDataOut into 3d data matrix (chan x time x epoch)
        sensorDataOut = reshape(sensorDataOut, sz(1), sz(2), []);
        
    otherwise
        error('Unknown montage type %s', montageType);
end

return