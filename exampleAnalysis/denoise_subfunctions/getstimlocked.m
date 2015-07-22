function sl = getstimlocked(data,which_freq)

% get the response amplitude at the stimulus locked frequency
% INPUT
% data   : raw time series [channels x time x epochs]
% freq   : index for stimulus locked frequency (either a number or a
%          struct)
% OUTPUT
% sl     : stimulus locked time series [epochs x channels]

% check input 
if notDefined('which_freq')
    which_freq = 12;
end

% Fourier transform data, take power of relevant frequency
spec = fft(data,[],2);
sl   = abs(squeeze(spec(:,which_freq,:)))' / size(data,2)*2;

return