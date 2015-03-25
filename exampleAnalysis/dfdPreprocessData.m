function sensorData = dfdPreprocessData(sensorDataIn, threshold)

%% Preprocess sensorData by identifying so-called 'bad' epochs and channels.
% These epochs and channels contain data whos variance is outside 
% multiples of the grand variance defined by the threshold.
%
% Example:
%  sensorData = dfdPreprocessData(sensorDataIn, [.01 10])

if notDefined('threshold'), threshold = [.05 20]; end

% This identifies any epochs whos variance is outside some multiple of the
% grand variance
bad_epochs = meg_find_bad_epochs(sensorDataIn, threshold);

% any epoch in which more than 10% of channels were bad should be removed
% entirely
epochs_to_remove = mean(bad_epochs,2)>.1;

% once we remove 'epochs_to_remove', check whether any channels have more
% than 10% bad epochs, and we will remove these
channels_to_remove = mean(bad_epochs(~epochs_to_remove,:),1)>.1;

bad_epochs(epochs_to_remove,:) = 1;
bad_epochs(:,channels_to_remove) = 1;

figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number');
fprintf('(dfdPreprocessData): Percentage of epochs removed is %5.2f \n', sum(sum(bad_epochs))/(size(sensorDataIn,2)*size(sensorDataIn,3))*100);

sensorData = meg_remove_bad_epochs(bad_epochs, sensorDataIn);

return

function bad_epochs = meg_find_bad_epochs(ts, thresh)
if ~exist('thresh', 'var') || isempty(thresh),
    thresh = [0.1 10];
end

var_matrix       = squeeze(nanvar(ts,[],1)); % var_matrix will be epochs x channels
var_grand_median = nanmedian(var_matrix(:)); % grand median

bad_epochs = var_matrix < thresh(1) * var_grand_median | ...
    var_matrix > thresh(2) * var_grand_median;
return


function ts = meg_remove_bad_epochs(bad_epochs, ts)
% epochs x channel
num_time_points = size(ts,1);
num_epochs      = size(ts,2);
num_channels    = size(ts,3);

ts = reshape(ts, [num_time_points, num_epochs*num_channels]);

ts(:, logical(bad_epochs(:))) = NaN;

ts = reshape(ts, [num_time_points, num_epochs, num_channels]);

return


