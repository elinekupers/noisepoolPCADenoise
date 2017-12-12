function [xFiltered, hipassfilter, freqs] = filterdata(x,fs,lcutoff,ftsd,excludef)
%
% Inputs
%   x:       3D array of data, [channels x time x epochs]
%   fs:      sampling frequency in Hz
%   lcutoff: low frequency cutoff for raised cosine high pass filter
%   ftsd:    width of the raised cosine high pass filter (in Hz)
%   excludef: band-stop frequencies to remove
%
% Outputs
%   xFiltered:      filtered data in same size and dimensions as x
%   hipassfilter:   amplitude spectrum of the filter
%   freqs:          frequencies (positive only)
%
% Example 1
%   x = randn(1, 1000,1);
%   [xFiltered, hipassfilter, freqs] = filterdata(x,1000);
%   figure; plot(abs(fft(x))); hold on, plot(abs(fft(xFiltered)), 'r'); xlim([0 200])
%
% Example 2
%   x = randn(1, 1000,1);
%   [xFiltered, hipassfilter, freqs] = filterdata(x,1000, [], [], 20:20:1000);
%   figure; plot(abs(fft(x))); hold on, plot(abs(fft(xFiltered)), 'r'); xlim([0 200])
%   figure; plot(abs(fft(xFiltered))./abs(fft(x))); xlim([0 200])

if notDefined('lcutoff'), lcutoff = 60; end
if notDefined('ftsd'),    ftsd = lcutoff; end

% get dimensions of data 
[nchan, xlen, nepochs] = size(x);
% create frequency axis
inds = 1:xlen; % inds = 1:floor(xlen/2)+1;
freqs= fs*linspace(0,1,length(inds)+1); 
inds = inds(1:round(end/2));
freqs = freqs(1:length(inds));

% create filter
hipassfilter = ones(1,xlen);
hipassfilter(freqs<lcutoff) = 0;
% create a smoothedge 
[xtbl,ytbl] = rcosFn(ftsd,lcutoff,[0,1]);
hipassfilter(1:length(freqs)) = interp1(xtbl,ytbl,freqs,'linear','extrap');
% also exclude certain frequencies 
if exist('excludef', 'var') && ~isempty(excludef)
    excludevec = ones(1,length(freqs));
    [~, ex_i]  = intersect(round(freqs), round(excludef));
    excludevec(ex_i) = 0;
    hipassfilter(1:length(freqs)) = hipassfilter(1:length(freqs)).*excludevec;
end
% mirror the smooth edge on the negative frequency side
if isodd(xlen)
    hipassfilter(xlen:-1:(round(xlen/2)+1)) = hipassfilter(2:round(xlen/2));
else
    hipassfilter(xlen:-1:(xlen/2+2)) = hipassfilter(2:xlen/2);
end
% make filter same dimensions as data 
hipassfilter2 = repmat(hipassfilter,[nchan,1,nepochs]);
% filter and transform back
xFiltered = real(ifft(fft(x,[],2).*hipassfilter2, [], 2));