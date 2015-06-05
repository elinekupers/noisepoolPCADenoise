function sl = getstimlocked(data,freq)
% get the stimulus locked frequency
% INPUT
% data   : raw time series [channels x time x epochs]
% freq   : index for stimulus locked frequency (either a number or a
%          struct)
% OUTPUT
% sl     : stimulus locked time series [epochs x channels]

% check input 
if isnumeric(freq)
    f = freq;
elseif isstruct(freq) && isfield(freq,'sl_i')
    f = freq.sl_i;
else
    error('input error: freq not recognized');
end

spec = fft(data,[],2);
sl   = abs(squeeze(spec(:,f,:)))' / size(data,2)*2;