function [pcnum,pcchan2,xvaltrend] = getpcchan(evalout,noisepool,pcn,pcstop)

% npcs x channels, averaged across non-noise channels
r2 = cat(1,evalout.r2);
% max cross validation for each channel
maxr2 = max(r2,[],1);
% exclude those in noisepool
maxr2(noisepool) = -inf;
% sort these
[~, idx] = sort(maxr2,'descend');
% pick the top x, determined by opt.pcn
pcchan = false(size(noisepool));
pcchan(idx(1:min(pcn,length(idx)))) = 1;
% take the average of these
xvaltrend = mean(r2(:,pcchan,1),2);
%xvaltrend = mean(r2(:,~noisepool),2);
% chosen the stopping point as within 95% of maximum
chosen  = choosepc(xvaltrend,pcstop);
pcchan2 = pcchan;
pcnum   = chosen;
