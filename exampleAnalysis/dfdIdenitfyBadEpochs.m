function [sensorData, okEpochs] = dfdIdenitfyBadEpochs(sensorDataIn, data_channels, threshold)
%


%%%%%%%%%%%%%%%%%%% THIS IS SSMEG PREPROC PART %%%%%%%%%%%%%%%%%%%
% Remove epochs with when the number of bad channels in that epoch exceeds threshold

% This identifies any epochs whos variance is outside some multiple of the
% grand variance
bad_epochs = meg_find_bad_epochs(sensorDataIn(:,:,data_channels), [.05 20]);

% any epoch in which more than 10% of channels were bad should be removed
% entirely
epochs_to_remove = mean(bad_epochs,2)>.1;

% once we remove 'epochs_to_remove', check whether any channels have more
% than 10% bad epochs, and we will remove these
channels_to_remove = mean(bad_epochs(~epochs_to_remove,:),1)>.1;

bad_epochs(epochs_to_remove,:) = 1;
bad_epochs(:,channels_to_remove) = 1;

figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number')

sensorData = meg_remove_bad_epochs(bad_epochs, sensorDataIn);

%%%%%%%%%%%%%%%%%%% THIS IS HELENA'S PART %%%%%%%%%%%%%%%%%%%
% sensor data is channels x time x epoch
badData = squeeze(isnan(sensorData(1,:,:)));
% average across channels
okEpochs = mean(badData)' < threshold;


% % remove bad epochs
% if opt.remove_badepochs
%     okEpochs = okEpochs & megIdenitfyBadEpochs(sensorData,0.5);
%     sensorData = sensorData(:,:,okEpochs);
%     epoch_conditions = epoch_conditions(okEpochs,:);
% end
% 
% % find bad channels - those where at least 50% are nan's 
% badChannels = megIdenitfyBadChannels(sensorData, 0.5);
% badChannels(98) = 1; % for now we add this in manually
% if opt.verbose
%     fprintf('\tbadChannels : '); 
%     fprintf('%g ', find(badChannels)');
%     fprintf('\n');
% end
% 
% % remove bad channels and replace remaining bad epochs 
% if opt.replace_badepochs
%     net=load('meg160xyz.mat');
%     sensorData = megReplaceBadEpochs(sensorData,net,[],opt.badepoch_avgchannum);
% end
% sensorData = sensorData(~badChannels,:,:);
% 
% 
% return
% 
% 
% function bad_epochs = meg_find_bad_epochs(ts, thresh)
% %Identify problem epochs in MEG time series. 
% %
% 
% if ~exist('thresh', 'var') || isempty(thresh), 
%     thresh = [0.1 10];
% end
% 
% var_matrix       = squeeze(nanvar(ts,[],1)); % var_matrix will be epochs x channels
% var_grand_median = nanmedian(var_matrix(:)); % grand median
% 
% bad_epochs = var_matrix < thresh(1) * var_grand_median | ...
%     var_matrix > thresh(2) * var_grand_median;

return

function ts = meg_remove_bad_epochs(bad_epochs, ts)
% epochs x channel
num_time_points = size(ts,1);
num_epochs      = size(ts, 2);
num_channels    = size(ts, 3);

ts = reshape(ts, [num_time_points, num_epochs*num_channels]);

ts(:, logical(bad_epochs(:))) = NaN;

ts = reshape(ts, [num_time_points, num_epochs, num_channels]);

return