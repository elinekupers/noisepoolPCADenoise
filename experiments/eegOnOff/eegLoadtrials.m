function [trialData,trialNums] = eegLoadtrials(pth, files, trials, channels, ...
    coherentAcrossEpochs, T)
% trialData = eegLoadtrials(pth, files, trials, channels, ...
%    coherentAcrossEpochs)
%
% Inputs
%
%
% Outputs
%
%
% Example

% Extract trial info from a sample file
hdr = eegGetTrialHeader(pth, files{1}(1));

% if isnumeric(trials),    % do nothing
% elseif strcmpi(trials, 'all'),    trials = 1:ntrials;
% else
%     warning('trials not specified. using all trials') %#ok<WNTAG>
%     trials = 1:ntrials;
% end

% How many channels?
if isnumeric(channels), % do nothing
elseif strcmpi(channels, 'all'), channels = 1:hdr.NmbChanEEG;
else
    warning('channels not specified. using all channels') %#ok<WNTAG>
    channels = 1:hdr.NmbChanEEG;
end

% How many epochs per trial?
if coherentAcrossEpochs, nmbEpochsPerTrial = 1;
else                     nmbEpochsPerTrial = (hdr.NmbEpochs-2)/T; end  % will this always be true? assumes epochs are 1 s long!

trialData = cell(1, numel(files));
trialNums = cell(1, numel(files));
% oneCondition = zeros(...
%     hdr.NmbChanEEG, ...                     number of channels
%     hdr.NmbSamples*2, ...                   number of time points per epoch
%     length(trials)*nmbEpochsPerTrial, ...   number of epochs
%     'double');                              % data format

% loop over conditions
for cond = 1:numel(files)
    
    % trialData{cond} = oneCondition;
            
    % loop over trials and extract data
    for trial = 1:length(files{cond})
        
        d = load(fullfile(pth, files{cond}(trial).name));
        
        % extract the signal matrix
        signal = d.RawTrial(:,1:hdr.NmbChanEEG);
        
        % recast as double precision
        signal = double(signal);

        % scale so that units are in units d.DataUnitStr        
        signal = bsxfun(@plus, signal, d.Shift(channels)');
        signal = bsxfun(@times, signal, d.Ampl(channels)');

        % reshape into time points x epochs x channels
        signal = reshape(signal, [], d.NmbEpochs, length(channels));
                
        % Epochs that are not ok get replaced with NaN
        for epoch = 1:d.NmbEpochs
            signal(:, epoch, ~d.IsEpochOK(epoch,:)) = NaN;
        end
         
        % Exclude 1st and last epoch of each trial
        signal = signal(:, 2:end-1, :);
                
        % Reshape into epochs of length T. If T is 1, this has no effect.
        signal = reshape(signal, hdr.NmbSamples*T, [], hdr.NmbChanEEG);
        
        % Scale into microVolts
        if strcmpi(hdr.DataUnitStr, 'volts'), signal = signal * 10E6; end
        
        % if we want the coherent spectrum...
        if coherentAcrossEpochs, signal = nanmean(signal,2); end
        
        % Permute signal into chan x time points x epoch
        signal = permute(signal, [3 1 2]);
        
        % tmp = cat(2, tmp, signal);
        theserows = (trial-1)*nmbEpochsPerTrial + (1:nmbEpochsPerTrial);
        trialData{cond}(:, :, theserows)  = signal; % tmp;
        trialNums{cond}(theserows,:) = trial*ones(size(theserows))';
    end
    
    
    
end

return