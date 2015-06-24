function [newdata, t] = hpf(data,sl_frequency)

if ~exist('sl_frequency', 'var'), sl_frequency = 12; end
ln = 60; % line noise

% harmonics of of the stimulus locked
tmp = [(1:40)*sl_frequency (1:10) * ln]; 
drop_frequencies  = sort(unique([tmp-1 tmp tmp+1]));
% high pass filter with cutoff of 62 Hz, sharp cutoff, and excluding
% harmonics
newdata = filterdata(data,1000,62,1,drop_frequencies);

t = size(newdata,2);