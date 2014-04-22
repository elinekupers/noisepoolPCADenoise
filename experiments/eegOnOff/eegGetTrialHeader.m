function hdr = eegGetTrialHeader(pth, file, T)

% T is the period for one epoch, in seconds. PowerDiva encodes a 12-second
% trial as 12 one-second epochs. We might choose to lump epochs into longer
% bigger bins, say 2-s, to get better frequency resolutiom. Or we may not.
if notDefined('T'), T = 1; end

% load a sample trial to extract some parameters that should be the same
% across trials
hdr = load(fullfile(pth, file.name));

% hdr.NmbEpochs;                          % number of epochs within a trial
% hdr.NmbChanEEG;                         % number of recording electrodes
% hdr.FreqHz;                             % sampling frequency
hdr.len        = size(hdr.RawTrial,1);    % num samples per trial
hdr.NmbSamples = hdr.len/hdr.NmbEpochs;   % number of samples per epoch
hdr.t = (1:hdr.NmbSamples*T)/hdr.FreqHz;  % time of one epoch
hdr.f = (0:length(hdr.t)-1)/max(hdr.t);   % frequencies for one epoch

return