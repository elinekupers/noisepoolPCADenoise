function ab = getbroadband(data,keep_frequencies,fs)

% get the broadband response
% INPUT
% data             : raw time series [channels x time x epochs]
% keep_frequencies : function handle to indicate which frequencies to use
%                       to calculate broadband
% fs               : sample rate (Hz)
%
% OUTPUT
% ab               : asynchronous broadband time series [epochs x channels]

% Define number of time points per epoch
n_time_points = size(data,2);

% Define frequencies to compute the broadband power
f = (0:n_time_points-1) * fs / n_time_points;
[~, f_inds] = intersect(f, keep_frequencies(f));

% Fourier transform the data, use only the relevant frequencies
spec = fft(data,[],2);
spec_amp = abs(spec(:,f_inds,:))/ size(data,2)*2; 

% take the mean across frequencies in log space 
ab = squeeze(exp(nanmean(log(spec_amp.^2),2)))';

return