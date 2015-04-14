function dfdCompariseDenoising(whichSubject, saveFigures, topChan, nTop)

if notDefined('whichSubject'); whichSubject = 1; end
if notDefined('saveFigures'); saveFigures = true; end
if notDefined('topchan'); topChan = true; end
if notDefined('nTop'); nTop = 10; end

% Load results with and without 3 channel denoising
tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_denoisedData.mat'),whichSubject)); 
noDenoise = tmp.results.origmodel;
megDenoise  = tmp.results.finalmodel;
a_noisepool = tmp.results.noisepool;

tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_denoisedData_w3chan_bb.mat'),whichSubject));
threechanDenoising = tmp.results.origmodel;
bothDenoise  = tmp.results.finalmodel;
b_noisepool = tmp.results.noisepool;


% Calculate SNR for all channels
%%
figure('position',[1,600,300,800]);
condNames    = {'Stim Full','Stim Left','Stim Right'};
condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
for icond = 1:3
    % get BB snr
    noDenoiseSNR(icond,:) = getsignalnoise(noDenoise,icond, 'SNR');
    megDenoiseSNR(icond,:) = getsignalnoise(megDenoise,icond, 'SNR');
    threechanDenoiseSNR(icond,:) = getsignalnoise(threechanDenoising,icond, 'SNR');
    bothDenoiseSNR(icond,:) = getsignalnoise(bothDenoise,icond, 'SNR');

end


%% If requested to only plot the top n channels, calculate those

if topChan
    %% Get top n channels from noisepool A
    % For noDenoise and MegDenoise first (they share the same noisepool)
    % max across 3 conditions
    topsnr(icond,:) = max([noDenoiseSNR(icond,:) megDenoiseSNR(icond,:)]);
    % max across before and after
    topsnr = max(topsnr);
    % exclude noise pool
    topsnr(a_noisepool) = -inf;
    % sort
    [~,idx] = sort(topsnr,'descend');
    % find the top n
    a_pcchan = false(size(a_noisepool));
    a_pcchan(idx(1:nTop))= 1;
    
    %% Get top n channels from noisepool B
    % For threechanDenoising and bothDenoise next (they share the same noisepool)
    % max across 3 conditions
    topsnr(icond,:) = max([threechanDenoiseSNR(icond,:) bothDenoiseSNR(icond,:)]);
    % max across before and after
    topsnr = max(topsnr);
    % exclude noise pool
    topsnr(b_noisepool) = -inf;
    % sort
    [~,idx] = sort(topsnr,'descend');
    % find the top n
    b_pcchan = false(size(b_noisepool));
    b_pcchan(idx(1:nTop))= 1;

end

% SNR No Denoising
subplot(4,1,1); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(noDenoiseSNR(nn,a_pcchan)),1),noDenoiseSNR(nn,a_pcchan),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(noDenoiseSNR(nn,a_pcchan),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR No Denoising', whichSubject));
    xlim([0.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
    ylabel('SNR')
end

% SNR Denoise only with the three environmental data channels
subplot(4,1,2); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(threechanDenoiseSNR(nn,b_pcchan))),threechanDenoiseSNR(nn,b_pcchan),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(threechanDenoiseSNR(nn,b_pcchan),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR 3 Channels', whichSubject));
    xlim([0.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
    ylabel('SNR')
end

% SNR Meg Denoising
subplot(4,1,3); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(megDenoiseSNR(nn,a_pcchan))),megDenoiseSNR(nn,a_pcchan),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(megDenoiseSNR(nn,a_pcchan),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR Meg Denoising', whichSubject));
    xlim([0.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
    ylabel('SNR')
end

% SNR of combined MEG Denoise and three environmental channels
subplot(4,1,4); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(bothDenoiseSNR(nn,b_pcchan))),bothDenoiseSNR(nn,b_pcchan),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(bothDenoiseSNR(nn,b_pcchan),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR Combined', whichSubject));
    xlim([.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
    ylabel('SNR')
end


if saveFigures
    if topChan
        figurewrite(fullfile(dfdRootPath,'figures',sprintf('figureCompareDenoisingTop%d',nTop)),[],0,'.',1);
    else
        figurewrite(fullfile(dfdRootPath,'figures','figureCompareDenoising'),[],0,'.',1);
    end
end
