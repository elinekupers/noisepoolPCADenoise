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
