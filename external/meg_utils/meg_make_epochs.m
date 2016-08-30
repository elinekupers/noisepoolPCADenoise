function [ts, conditions] = meg_make_epochs(raw_ts, trigger, epoch_time, fs, which_data)
% Slice time series matrix (samples x channels) into 3D epoched array
% (samples x epoch x channel) based on trigger times.
%
%[ts, conditions] = meg_make_epochs(raw_ts, trigger, epoch_time, fs)
%
% Inputs
%   raw_ts:       time series matrix (samples x channels)
%   trigger:      a vector of stimulus onsets equal in length to first
%                   dimension of raw_ts. Should be all zeros except an
%                   integer to indicate trial onset. These integer values
%                   correspond to the condition number.
%   epoch_times:  a 2 vector of start and end time of the epochs 
%                   (in seconds) relative to the trial onset
%   fs:           sampling rate (Hz)
%   which_data:   string to define which data is used, if 'eye',
%                   then triggers are already defined as timepoints   
%
% Outputs
%   ts:           3D array containing epoched time series (samples by
%                   epoch x channel)
%   conditions:   vector equal in length to the number of epochs. each
%                   entry is the condition number (trigger value) for that
%                   epoch

%% Parameters

if ~exist('which_data','var') 
    onset_times = find(trigger);
    which_data = 'meg';
elseif which_data == 'eye';
    onset_times = trigger;
    if onset_times(1) == 0; onset_times(1) = 1; end % First trigger can't be 0.
    last_idx = length(raw_ts);    
end

epoch_samples = round(epoch_time * fs); %epoch length in samples

inds = bsxfun(@plus,onset_times,(epoch_samples(1):epoch_samples(2)));

% if the onset times + number of timepoints are longer than the actual raw
% ts, we have to discard those onset times
if which_data == 'eye';
    if ~isempty(find(inds(:)>last_idx))
        discarded_onsets = ceil(length(find(inds(:)>last_idx))/diff(epoch_samples));
        inds = inds(1:end-discarded_onsets,:);
    end
end

ts   = raw_ts(inds,:); 
ts   = reshape(ts,size(inds,1),size(inds,2),size(raw_ts,2));
ts   = permute(ts, [2 1 3]);

if strcmp(which_data,'meg'); conditions = trigger(onset_times); 
else conditions = []; end
return
