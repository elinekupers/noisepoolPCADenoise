function newdata = bbFilter(data,keep_frequencies, fs)

% [newdata, t] = bbFilter(data,keep_frequencies, fs)
%
% Function to filter data by removing all frequencies except those that are
% used to compute broadband power
%
% INPUTS
% data              : Epoched timeseries [channels x timepoints x epochs]
% keep_frequencies  : A vector of frequencies to keep
% fs                : sampling frequency [default = 1000]

% OUTPUTS
% newdata       :   Filtered data [channels x timepoints x epochs]            

% Check inputs
if ~exist('fs', 'var') || isempty(fs), fs = 1000; end

freq = 0:fs/size(data,2):fs/2;
drop_f = setdiff(freq, keep_frequencies);

if isempty(intersect(freq, keep_frequencies)), 
    error('No frequencies found to keep. All frequencies will be removed by filter.')
end

newdata = filterdata(data,fs,1,1,drop_f);

return
