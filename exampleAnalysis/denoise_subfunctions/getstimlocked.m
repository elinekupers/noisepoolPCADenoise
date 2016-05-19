function sl = getstimlocked(data,indexed_frequency)

% get the response amplitude at the stimulus locked frequency
% INPUT
% data   : raw time series [channels x time x epochs]
% freq   : index for stimulus locked frequency (either a number or a
%          struct)
% OUTPUT
% sl     : stimulus locked time series [epochs x channels]

% check input 
if notDefined('which_freq')
    indexed_frequency = 13; % corresponds to 12 Hz at 1 Hz spacing
end

% Fourier transform data, take amplitude of relevant frequency
spec = fft(data,[],2);
sl   = abs(squeeze(spec(:,indexed_frequency,:)))' / size(data,2)*2;

return