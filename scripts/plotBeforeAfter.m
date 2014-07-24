function plotBeforeAfter(allresults, ...
    whichfun, doTop10, plotType, whichConds, funXchan)
% function to plot denoising results before and after denoising, across
% subjects
% allresults is a cell of length nsessions, each containing results struct
% from denoisedata
% rest of the arguments are optional 

if notDefined('whichfun'), whichfun = 1;    end  % which function, if results contains more than one
if notDefined('doTop10'),  doTop10  = true; end  % top 10 or all non-noise
if notDefined('plotType'), plotType = 'snr'; end % snr, s, or n 
if notDefined('whichConds'), whichConds = 1:3; end % 1:3, or 1, 2, 3 (full, right, left)
if notDefined('funXchan'),  funXchan = @mean; end  % how to aggregate across chan (e.g.@mean or @median)

nsess = length(allresults);

% values before and after denoising, nsess x 1 
vals1all = [];
vals2all = [];

for k = 1:nsess
    results = allresults{k};
    
    % toggle between the top10 channels or all non-noise channels
    if doTop10
        try
            pcchan = results.pcchan{whichfun};
        catch exception
            disp(exception.message);
            fprintf('\tsetting pcchan to top10 of final model\n');
            finalsnr = getsignalnoise(results.finalmodel(whichfun));
            finalsnr(results.noisepool) = -inf;
            [~,idx] = sort(finalsnr,'descend');
            pcchan = false(size(results.noisepool));
            pcchan(idx(1:10))= 1;
        end
    else
        pcchan = ~results.noisepool; 
    end
    
    % nconds x nchannels
    ab_signal1 = abs(results.origmodel(whichfun).beta_md(:,pcchan));
    ab_noise1  = results.origmodel(whichfun).beta_se(:,pcchan);
    ab_signal2 = abs(results.finalmodel(whichfun).beta_md(:,pcchan));
    ab_noise2  = results.finalmodel(whichfun).beta_se(:,pcchan);
    ab_snr1    = ab_signal1./ab_noise1;
    ab_snr2    = ab_signal2./ab_noise2;
    
    % choose the values we want 
    switch plotType
        case 'snr'
            vals1 = max(ab_snr1(whichConds,:),[],1);
            vals2 = max(ab_snr2(whichConds,:),[],1);
        case 's'
            vals1 = max(ab_signal1(whichConds,:),[],1);
            vals2 = max(ab_signal2(whichConds,:),[],1);
        case 'n'
            vals1 = mean(ab_noise1(whichConds,:),1);
            vals2 = mean(ab_noise2(whichConds,:),1);
    end
    vals1all = cat(1, vals1all, funXchan(vals1));
    vals2all = cat(1, vals2all, funXchan(vals2));
end

% plot across subjects/sessions 
colors = copper(nsess);
hold on;
for k = 1:nsess
    plot(1:2, [vals1all(k),vals2all(k)], 'o-', 'color', colors(k,:), 'linewidth',2);
end
xlim([0,3]);
set(gca,'xtick',1:2,'xticklabel',{'Before','After'}); 
ylabel(upper(plotType));
makeprettyaxes(gca,14,14); axis square;
% make title 
if length(whichConds) > 1, condName = 'all'; else condName = num2str(whichConds); end
if doTop10, st = 'Top 10'; else st = 'Non Noise'; end
ttl = sprintf('Cond: %s, %s across %s channels', condName, func2str(funXchan), st);
title(ttl);