function dfdCompariseDenoising(whichSubject, saveFigures)

if notDefined('whichSubject'); whichSubject = 1; end
if notDefined('saveFigures'); saveFigures = false; end

% Load results with and without 3 channel denoising
tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_denoisedData.mat'),whichSubject)); 
noDenoise = tmp.results.origmodel;
megDenoise  = tmp.results.finalmodel;

tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_denoisedData_w3chan_bb.mat'),whichSubject));
threechanDenoising = tmp.results.origmodel;
bothDenoise  = tmp.results.finalmodel;


% Calculate SNR for all channels
%%
figure('position',[1,600,300,800]);
condNames    = {'Stim Full','Stim Left','Stim Right'};
condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
for icond = 1:3
    % get BB snr
    noDenoiseSNR(icond,:) = getsignalnoise(noDenoise,icond, 'SNR');
    megDenoiseSNR(icond,:) = getsignalnoise(megDenoise,icond, 'SNR');
    threechanDenoisingSNR(icond,:) = getsignalnoise(threechanDenoising,icond, 'SNR');
    bothDenoiseSNR(icond,:) = getsignalnoise(bothDenoise,icond, 'SNR');
end


% SNR No Denoising
subplot(4,1,1); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(noDenoiseSNR(nn,:)),1),noDenoiseSNR(nn,:)','o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(noDenoiseSNR(nn,:),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR No Denoising', whichSubject));
    xlim([0.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
end

% SNR Denoise only with the three environmental data channels
subplot(4,1,2); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(threechanDenoisingSNR(nn,:))),threechanDenoisingSNR(nn,:),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(threechanDenoisingSNR(nn,:),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR 3 Channels', whichSubject));
    xlim([0.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
end

% SNR Meg Denoising
subplot(4,1,3); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(megDenoiseSNR(nn,:))),megDenoiseSNR(nn,:),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(megDenoiseSNR(nn,:),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR Meg Denoising', whichSubject));
    xlim([0.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
end

% SNR of combined MEG Denoise and three environmental channels
subplot(4,1,4); cla; hold on;
for nn = 1:3
    plot(nn*ones(length(bothDenoiseSNR(nn,:))),bothDenoiseSNR(nn,:),'o','color',condColors(nn,:));
    plot(nn*ones(1,1),mean(bothDenoiseSNR(nn,:),2),'kx');
    axis tight;
    title(sprintf('S%d : SNR Combined', whichSubject));
    xlim([.5,3.5]); ylim([0,15]);
    makeprettyaxes(gca,9,9);
end


if saveFigures
    figurewrite(fullfile(dfdRootPath,'figures','figureCompareDenoising'),[],0,'.',1);
end
