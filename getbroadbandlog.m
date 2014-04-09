function ab = getbroadbandlog(data,freq)

% get the stimulus locked frequency
% INPUT
% data   : raw time series [channels x time x epochs]
% freq   : indicies for frequencies to compute (either a vector or a
%          struct)
% OUTPUT
% sl     : broadband time series [epochs x channels]

% check input 
if isnumeric(freq)
    f = freq;
elseif isstruct(freq) && isfield(freq,'sl_i')
    f = freq.ab_i;
else
    error('input error: freq not recognized');
end

spec = fft(data,[],2);
spec_amp = abs(spec(:,f,:))/ size(data,2)*2;
% this makes sure we're not taking log of zero
spec_amp(spec_amp==0) = eps;
% take the mean across frequencies in log space 
ab = squeeze(nanmean(log(spec_amp.^2),2))';

