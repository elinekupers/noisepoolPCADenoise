function allvals = getSNRgrid(allresults,npools,npcs,...
    whichfun, doTop10, plotType, whichConds, funXchan, plotFlg)
% see also megDenoiseGridScript for usage

if notDefined('whichfun'), whichfun = 1;    end  % which evalfun, if results contains more than one
if notDefined('doTop10'),  doTop10  = true; end  % top 10 or all non-noise
if notDefined('plotType'), plotType = 'snr'; end % snr, s, or n
if notDefined('whichConds'), whichConds = 1:3; end % 1:3 (max of 3), or 1, 2, 3 (full, right, left)
if notDefined('funXchan'),  funXchan = @mean; end  % how to aggregate across chan (e.g.@mean or @median)
if notDefined('plotFlg'),  plotFlg = true; end  % plot?

allvals = nan(length(npools),length(npcs));
for np = 1:length(npools)
    for nc = 1:length(npcs)
        results = allresults(np,nc);
        if isempty(results.finalmodel), continue; end
        
        % toggle between the top10 channels or all non-noise channels
        if doTop10
            try
                pcchan = results.pcchan{whichfun};
            catch exception
                %disp(exception.message);
                %fprintf('\tsetting pcchan to top10 of final model\n');
%                 finalsnr = getsignalnoise(results.finalmodel(whichfun));
%                 finalsnr(results.noisepool) = -inf;
%                 [~,idx] = sort(finalsnr,'descend');
%                 pcchan = false(size(results.noisepool));
%                 pcchan(idx(1:10))= 1;
                pcchan = getTop10(results,whichfun);
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
                vals1 = 1./mean(ab_noise1(whichConds,:),1);
                vals2 = 1./mean(ab_noise2(whichConds,:),1);
            case 'r2'
                vals1 = results.origmodel(whichfun).r2(:,pcchan);
                vals2 = results.finalmodel(whichfun).r2(:,pcchan);
        end
        
        % take the median or mean across channels, and find the difference
        % between before and after 
        allvals(np,nc) = funXchan(vals2)-funXchan(vals1);
    end
end

if plotFlg
    plotSNRgrid(allvals,npools,npcs, whichfun, doTop10, plotType, whichConds, funXchan);
end