function sensorDataOut = eegReplaceBadEpochs(sensorDataIn, montageType)
% sensorDataOut = eegReplaceBadChannels(sensorDataIn, montageType)


% Replace channels with mean of neighboring channels
if notDefined('montageType'), montageType = 'egi'; end
if ~exist('defaultSphereNet.mat', 'file'), addpath('mrCurrent'); end

sensorDataOut = cell(size(sensorDataIn));

for ii = 1:numel(sensorDataIn)
    % sensor data is channels x time x epoch
    sz = size(sensorDataIn{ii});
    
    switch lower(montageType)
        case {'128', 'egi', 'hd'}
            % we expect 128 channels
            assert(sz(1) == 128);
            
            % load the xyz positions of the sensor
            net=load('defaultSphereNet.mat');
            
            % calculate a distance matrix
            distances = squareform(pdist(net.xyz), 'tomatrix');
            
            % for each channel, find neighbors within X mm
            nearest = distances < 4 & distances > 0;
                        
            % the sum of the columns should be one
            nearest = bsxfun(@rdivide, nearest, sum(nearest));
            
            % eliminate any NaNs (channels with no neighbors)
            nearest(isnan(nearest)) = 0;
            
            % reshape 3d sensor data into matrix (2d) in order to multiply
            sensorDataIn{ii} = reshape(sensorDataIn{ii},  sz(1), []);
            
            % Initialize out data
            sensorDataOut{ii} = sensorDataIn{ii};
            
            
            % Replace bad channels with the mean of the neighbors

            for jj = 1:128
                inds = find(nearest(:, jj));
                badPoints = isnan(sensorDataOut{ii}(jj,:));
                tmp = sensorDataIn{ii}(inds,badPoints);
                tmp(isnan(tmp)) = 0;
                sensorDataOut{ii}(jj,badPoints) = (tmp' * nearest(inds, jj))';
            end
            % reshape 2d sensorDataOut into 3d data matrix (chan x time x epoch)
            sensorDataOut{ii} = reshape(sensorDataOut{ii}, sz(1), sz(2), []);
            
        otherwise
            error('Unknown montage type %s', montageType);
    end
    
end

return