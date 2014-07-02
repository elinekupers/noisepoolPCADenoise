function [xFiltered, hipassfilter, freqs] = filterdata(x,fs,lcutoff,ftsd,excludef)
%fs = sampling frequency
%lcutoff = lower frequency cutoff
if notDefined('lcutoff'), lcutoff = 60; end
if notDefined('ftsd'),    ftsd = lcutoff; end

% get dimensions of data 
[nchan, xlen, nepochs] = size(x);
% create frequency axis
inds = 1:floor(xlen/2)+1;
freqs= fs/2*linspace(0,1,length(inds));
% create filter
hipassfilter = ones(1,xlen);
hipassfilter(freqs<lcutoff) = 0;
% create a smoothedge 
[xtbl,ytbl] = rcosFn(ftsd,lcutoff,[0,1]);
hipassfilter(1:length(freqs)) = interp1(xtbl,ytbl,freqs,'linear','extrap');
% also exclude certain frequencies 
if ~notDefined('excludef')
    excludevec = ones(1,length(freqs));
    [~, ex_i]  = intersect(round(freqs), excludef);
    excludevec(ex_i) = 0;
    hipassfilter(1:length(freqs)) = hipassfilter(1:length(freqs)).*excludevec;
end
% mirror the smoothedge on the negative frequency side
if isodd(xlen)
    hipassfilter(xlen:-1:(round(xlen/2)+1)) = hipassfilter(2:round(xlen/2));
else
    hipassfilter(xlen:-1:(xlen/2+2)) = hipassfilter(2:xlen/2);
end
% make filter same dimensions as data 
hipassfilter2 = repmat(hipassfilter,[nchan,1,nepochs]);
% filter and transform back
xFiltered = real(ifft(fft(x,[],2).*hipassfilter2, [], 2));