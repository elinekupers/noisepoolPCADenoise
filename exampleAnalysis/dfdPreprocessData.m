function sensorData = dfdPreprocessData(sensorDataIn, threshold)

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

figure; imagesc(bad_epochs); xlabel('channel number'); ylabel('epoch number')

sensorData = meg_remove_bad_epochs(bad_epochs, sensorDataIn);











%% Remove bad epochs



%% Save prepocessed data



