function trial_data = epoch2trial(trial_data, epoched_data)

num_channels = size(epoched_data,1);
num_trials = numel(trial_data);
epoched_data = reshape(epoched_data, num_channels,[], num_trials);

num_trial_time_points = size(epoched_data, 2);
for ii = 1:num_trials
   trial_data{ii}(:, 1:num_trial_time_points) = epoched_data(:, :, ii);
end

end
