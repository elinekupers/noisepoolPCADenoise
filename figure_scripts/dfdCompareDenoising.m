function dfdCompareDenoising(whichSubjects, saveFigures, topChan, nTop)

if notDefined('whichSubject'); whichSubjects = 1:8; end
if notDefined('saveFigures'); saveFigures = true; end
if notDefined('topchan'); topChan = true; end
if notDefined('nTop'); nTop = 10; end

allResults = [];

for whichSubject = whichSubjects
    % Load results with and without 3 channel denoising
    tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_denoisedData_bb.mat'),whichSubject)); 
    noDenoise = tmp.results.origmodel;
    megDenoise  = tmp.results.finalmodel;
    a_noisepool = tmp.results.noisepool;

    tmp = load(sprintf(fullfile(dfdRootPath, 'data', 's0%d_denoisedData_w3chan_bb.mat'),whichSubject));
    threechanDenoising = tmp.results.origmodel;
    bothDenoise  = tmp.results.finalmodel;
    b_noisepool = tmp.results.noisepool;


% Calculate SNR for all channels
%%
% figure('position',[1,600,300,800]);
    noDenoiseSNR = [];
    megDenoiseSNR = [];
    threechanDenoiseSNR = [];
    bothDenoiseSNR = [];
    
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
        topsnr = max([noDenoiseSNR; megDenoiseSNR]);
           
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
        topsnr = max([threechanDenoiseSNR; bothDenoiseSNR]);
        
        % exclude noise pool
        topsnr(b_noisepool) = -inf;
        % sort
        [~,idx] = sort(topsnr,'descend');
        % find the top n
        b_pcchan = false(size(b_noisepool));
        b_pcchan(idx(1:nTop))= 1;


    end

    % Check whether top n channels of two loaded files are the same
    assert(isequal(a_pcchan, b_pcchan))
    pcchan = a_pcchan;
    
    allResults = {noDenoiseSNR,megDenoiseSNR,threechanDenoiseSNR,bothDenoiseSNR,pcchan};
    
    % Put everything in one array (conditions by 5 times SNR beta's)
    results_null{whichSubject} = catcell(1,allResults);
    
    
end

%% Compute difference in pre-post SNR

    
    % Get the data for each kind of condition, for each subject so we can
    % take the mean/median across subjects
%     allsubjects_noDenoiseSNR

snr_pre1  = []; % No denoise
snr_post1 = []; % 3 Channel environmental denoise

snr_diff1 = []; % No denoise vs 3 Channel environmental denoise

snr_pre2  = []; % No denoise
snr_post2 = []; % MEG Denoise

snr_diff2 = []; % No denoise vs MEG Denoise

snr_pre3  = []; % No denoise
snr_post3 = []; % Combination of MEG Denoise & 3 Channel Environmental d.

snr_diff3 = []; % No denoise vs Combination of MEG Denoise & 3 Channel Environmental d.


% only full for now
for whichSubject = whichSubjects
    for icond = 1:3
        pcchan = find(results_null{whichSubject}(end,:));
        
        % Environmental denoising
        snr_pre1(whichSubject,icond,:) = results_null{whichSubject}(icond,pcchan);
        snr_post1(whichSubject,icond,:) = results_null{whichSubject}(6+icond,pcchan);
        
        
        % MEG denoising
        snr_pre2(whichSubject,icond,:) = results_null{whichSubject}(icond,pcchan);
        snr_post2(whichSubject,icond,:) = results_null{whichSubject}(3+icond,pcchan);
        
        
        % BOTH denoising        
        snr_pre3(whichSubject,icond,:) = results_null{whichSubject}(icond,pcchan);
        snr_post3(whichSubject,icond,:) = results_null{whichSubject}(9+icond,pcchan);
        
        
    end
    
    
end

    % subjects by conditions (8*3)
    snr_diff1 = mean(snr_post1 - snr_pre1,3); %mean across diff of pre vs post for top 10 channels
    snr_diff2= mean(snr_post2 - snr_pre2,3);
    snr_diff3 = mean(snr_post3 - snr_pre3,3);
    
% For one subject only / obsolete now
%     % No denoise vs 3 Channel environmental denoise
%     snr_pre1  = noDenoiseSNR(icond,pcchan);
%     snr_post1 = threechanDenoiseSNR(icond,pcchan);
%     
%     snr_diff(whichSubject,1,icond,:) = snr_post1 - snr_pre1;
%     
%     % No denoise vs MEG Denoise
%     snr_pre2  = noDenoiseSNR(icond,pcchan);
%     snr_post2 = megDenoiseSNR(icond,pcchan);
%     
%     snr_diff(whichSubject,2,icond,:) = snr_post2 - snr_pre2;
%     
%     % No denoise vs Combination of MEG Denoise & 3 Channel Environmental d.
%     snr_pre3  = noDenoiseSNR(icond,pcchan);
%     snr_post3 = bothDenoiseSNR(icond,pcchan);
%     
%     snr_diff(whichSubject,3,icond,:) = snr_post3 - snr_pre3;

% subjects x conditions x 
snr_diff = cat(3,snr_diff1,snr_diff2,snr_diff3);

condNames    = {'Stim Full','Stim Left','Stim Right'};
condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;

fH = figure('position',[0,300,700,200]);
% define what the different conditions are 
types = {'No denoise vs environmental','No denoise vs MEG denoise','No Denoise vs Combined denoise'};
% % re-arrange the order of the bars 
% neworder = [1,2,3];
% newtypes = types(neworder);
% 
% snr_diff2 = snr_diff(:,neworder,:);
% nnull = length(types);
% fH = figure('position',[0,300,200,200]);
nnull = length(types);

for icond = 1:3
%     for ncontrol = 1:size(snr_diff,3)
        subplot(1,3,icond);
        % mean and sem across subjects
        mn  = mean(snr_diff(:,icond,:),1);
        sem = std(snr_diff(:,icond,:),[],1)/sqrt(length(whichSubjects));
        bar(1:nnull, squeeze(mn),'EdgeColor','none','facecolor',condColors(icond,:)); hold on
        errorbar2(1:nnull,squeeze(mn),squeeze(sem),1,'-','color',condColors(icond,:));
        % format figure and make things pretty
        set(gca,'xlim',[0.2,nnull+0.8],'ylim',[-1,5]);
        makeprettyaxes(gca,9,9);
        %     set(gca,'XLabel',types{1:3})
%     end
end



if saveFigures
    if topChan
        figurewrite(fullfile(dfdRootPath,'figures',sprintf('figureCompareDenoisingTop%d',nTop)),[],0,'.',1);
    else
        figurewrite(fullfile(dfdRootPath,'figures','figureCompareDenoising'),[],0,'.',1);
    end
end




% %% Plot all channels for the three conditions w/4 different analysis
% % SNR No Denoising
% subplot(4,1,1); cla; hold on;
% for nn = 1:3
%     plot(nn*ones(length(noDenoiseSNR(nn,a_pcchan)),1),noDenoiseSNR(nn,a_pcchan),'o','color',condColors(nn,:));
%     plot(nn*ones(1,1),mean(noDenoiseSNR(nn,a_pcchan),2),'kx');
%     axis tight;
%     title(sprintf('S%d : SNR No Denoising', whichSubject));
%     xlim([0.5,3.5]); ylim([0,15]);
%     makeprettyaxes(gca,9,9);
%     ylabel('SNR')
% end
% 
% % SNR Denoise only with the three environmental data channels
% subplot(4,1,2); cla; hold on;
% for nn = 1:3
%     plot(nn*ones(length(threechanDenoiseSNR(nn,b_pcchan))),threechanDenoiseSNR(nn,b_pcchan),'o','color',condColors(nn,:));
%     plot(nn*ones(1,1),mean(threechanDenoiseSNR(nn,b_pcchan),2),'kx');
%     axis tight;
%     title(sprintf('S%d : SNR 3 Channels', whichSubject));
%     xlim([0.5,3.5]); ylim([0,15]);
%     makeprettyaxes(gca,9,9);
%     ylabel('SNR')
% end
% 
% % SNR Meg Denoising
% subplot(4,1,3); cla; hold on;
% for nn = 1:3
%     plot(nn*ones(length(megDenoiseSNR(nn,a_pcchan))),megDenoiseSNR(nn,a_pcchan),'o','color',condColors(nn,:));
%     plot(nn*ones(1,1),mean(megDenoiseSNR(nn,a_pcchan),2),'kx');
%     axis tight;
%     title(sprintf('S%d : SNR Meg Denoising', whichSubject));
%     xlim([0.5,3.5]); ylim([0,15]);
%     makeprettyaxes(gca,9,9);
%     ylabel('SNR')
% end
% 
% % SNR of combined MEG Denoise and three environmental channels
% subplot(4,1,4); cla; hold on;
% for nn = 1:3
%     plot(nn*ones(length(bothDenoiseSNR(nn,b_pcchan))),bothDenoiseSNR(nn,b_pcchan),'o','color',condColors(nn,:));
%     plot(nn*ones(1,1),mean(bothDenoiseSNR(nn,b_pcchan),2),'kx');
%     axis tight;
%     title(sprintf('S%d : SNR Combined', whichSubject));
%     xlim([.5,3.5]); ylim([0,15]);
%     makeprettyaxes(gca,9,9);
%     ylabel('SNR')
% end
% 
% 
% if saveFigures
%     if topChan
%         figurewrite(fullfile(dfdRootPath,'figures',sprintf('figureCompareDenoisingTop%d',nTop)),[],0,'.',1);
%     else
%         figurewrite(fullfile(dfdRootPath,'figures','figureCompareDenoising'),[],0,'.',1);
%     end
% end
