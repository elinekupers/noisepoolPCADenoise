function [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorDataIn, varThreshold, ...
    badChannelThreshold, badEpochThreshold)
% Preprocess MEG data
%
%sensorData = dfdPreprocessData(sensorDataIn, varThreshold, ...
%    badChannelThreshold, badEpochThreshold)
%
% INPUTS
%   sensorDataIn: 3D array, time points x epochs x channels
%
%   varThreshold: Vector of length 2 ([min max]) to indicate variance
%                   threshold. For any channel in any give epoch, if the
%                   variance in the time series is outside [min max] * the
%                   median variance across all channels and all epochs,
%                   then label that channel during that epoch as 'bad'
%                       Default = [.05 20]
%
%  badChannelThreshold: Fraction ([0 1]). If more than this fraction of
%                   channels in any given epoch is labeled 'bad', then
%                   label all channels for this epoch as 'bad'.
%                       Default = 0.2
%
%  badEpochThreshold: Fraction ([0 1]). If more than this fraction of
%                   epochs for any given channel is labeled 'bad', then
%                   label all epochs for this  channel as 'bad'
%                       Default = 0.2
%
% OUTPUTS
%   sensorData:     Same as sensorDataIn (3D array, time points x epochs x
%                   channels), except that for which 'bad', data has been
%                   replaced with NaN
%
%
% Example:
%  sensorData = dfdPreprocessData(sensorDataIn, [.01 10], .2, .2);

if notDefined('varThreshold'), varThreshold = [.05 20]; end
if notDefined('badChannelThreshold'), badChannelThreshold = .2; end
if notDefined('badEpochThreshold'), badEpochThreshold = .2; end

% This identifies any epochs whos variance is outside some multiple of the
% grand variance
outliers = meg_find_bad_epochs(sensorDataIn, varThreshold);

% any epoch in which more than 10% of channels were bad should be removed
% entirely
badEpochs = mean(outliers,2)>badChannelThreshold;

% once we remove 'badEpochs', check whether any channels have more
% than 10% bad epochs, and we will remove these
badChannels = mean(outliers(~badEpochs,:),1)>badEpochThreshold;

outliers(badEpochs,:)   = 1;
outliers(:,badChannels) = 1;

figure; imagesc(outliers); 
xlabel('channel number'); ylabel('epoch number'); title('Bad channels / epochs')
fprintf('(dfdPreprocessData): %5.2f%% of epochs removed\n', sum(sum(outliers))/(size(sensorDataIn,2)*size(sensorDataIn,3))*100);

sensorData = dfdChannelRepair(sensorDataIn, outliers, 'nearest');


% sensorData = meg_remove_bad_epochs(outliers, sensorDataIn);
% 
% % Remove bad channels and bad epochs
% sensorData = sensorData(:, ~badEpochs, ~badChannels);
% 
% % Permute sensorData for denoise code, which expects channel x time x epochs
% sensorData = permute(sensorData, [3,1,2]); % channel x time samples x epoch

return

function outliers = meg_find_bad_epochs(ts, thresh)
if ~exist('thresh', 'var') || isempty(thresh),
    thresh = [0.1 10];
end

var_matrix       = squeeze(nanvar(ts,[],1)); % var_matrix will be epochs x channels
var_grand_median = nanmedian(var_matrix(:)); % grand median

outliers = var_matrix < thresh(1) * var_grand_median | ...
    var_matrix > thresh(2) * var_grand_median;
return


function ts = meg_remove_bad_epochs(outliers, ts)
% epochs x channel
num_time_points = size(ts,1);
num_epochs      = size(ts,2);
num_channels    = size(ts,3);

ts = reshape(ts, [num_time_points, num_epochs*num_channels]);

ts(:, logical(outliers(:))) = NaN;

ts = reshape(ts, [num_time_points, num_epochs, num_channels]);

return

% function ts = channelrepair(ts, 'method')

% Use as
%   [interp] = ft_channelrepair(cfg, data)
%
% The configuration must contain
%   cfg.method         = 'nearest', 'average', 'spline' or 'slap' (default='nearest')
%   cfg.badchannel     = cell-array, see FT_CHANNELSELECTION for details
%   cfg.missingchannel = cell-array, see FT_CHANNELSELECTION for details
%   cfg.neighbours     = neighbourhoodstructure, see also FT_PREPARE_NEIGHBOURS
%   cfg.trials         = 'all' or a selection given as a 1xN vector (default = 'all')
%   cfg.lambda         = regularisation parameter (default = 1e-5, not for method 'distance')
%   cfg.order          = order of the polynomial interpolation (default = 4, not for method 'distance')

