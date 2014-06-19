function [pcnum,pcchan,xvaltrend] = getpcchan(evalout,noisepool,pcn,pcstop,pcselmethod)

if notDefined('pcselmethod'), pcselmethod = 'r2'; end
% npcs2try x channels, averaged across non-noise channels
switch pcselmethod
    case 'r2'
        metric = cat(1,evalout.r2);
    case 'snr'
        metric = max(abs(cat(3,evalout.beta_md)),[],1) ./ mean(cat(3,evalout.beta_se),1);
        metric = squeeze(metric)';
end
% max value across npcs2try for each channel
maxmetric = max(metric,[],1);
% exclude those in noisepool
maxmetric(noisepool) = -inf;
% sort these
[~, idx] = sort(maxmetric,'descend');
% pick the top x, determined by opt.pcn
pcchan = false(size(noisepool));
pcchan(idx(1:min(pcn,length(idx)))) = 1;
% take the average of these
xvaltrend = mean(metric(:,pcchan),2);
%xvaltrend = mean(r2(:,~noisepool),2);
% chosen the stopping point as within 95% of maximum
chosen  = choosepc(xvaltrend,pcstop);
pcnum   = chosen;

