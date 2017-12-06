function fH = plotEpochLengthVersusNPCsVersusSNR(whichSubjects,dataAll,epochDurs,npcs,saveFigures,figureDir)

snr_diff = zeros(length(whichSubjects),length(epochDurs),length(npcs),3);

for whichSubject = whichSubjects
    
    % get all the results for a particular number of pc's removed
    for jj = 1:length(npcs)
        
        results_all = catcell(1,dataAll{whichSubject}{1}.allResults(:,jj));
        
        % get top 10 channels
%         pcchan = getTop10(results_all(1));
        
        % compute the difference between pre and post
        for nn = 1:length(epochDurs)
            pcchan = getTop10(results_all(nn));
            
            bb_signal_pre = results_all(nn).origmodel(1).beta_mn(:,pcchan);
            bb_noise_pre  = results_all(nn).origmodel(1).beta_sd(:,pcchan);
            
            bb_signal_post = results_all(nn).finalmodel(1).beta_mn(:,pcchan);
            bb_noise_post  = results_all(nn).finalmodel(1).beta_sd(:,pcchan);
            
            bb_snr_pre     = bb_signal_pre./bb_noise_pre;
            bb_snr_post    = bb_signal_post./bb_noise_post;
            
            
            snr_diff(whichSubject,nn,jj,:) = mean((bb_snr_post-bb_snr_pre),2);
            
        end
    end
end

%% Plot
fH = figure('position',[0,300,300,600]);
epochDurs = [1,3,6,12,24,36,72,1080];
clims = [[0,4];[0,3];[0,3]];
% plot each condition as a separate panel
for icond = 1:3
    subplot(3,1,icond);
    % plot the mean across subject as one square
    sessionAvg = squeeze(mean(snr_diff(:,:,:,icond),1));
    imagesc(sessionAvg',clims(icond,:));
    % set up axes and format figure
    set(gca,'ydir','normal');
    set(gca,'xtick',1:length(epochDurs),'ytick',1:length(npcs),...
        'xticklabel',cellstr(num2str(epochDurs','%d')),'yticklabel',cellstr(num2str(npcs','%d')));
    axis image;
    set(gca,'tickdir','out'); box off
    ch = colorbar; makeprettyaxes(ch,9,9);
    xlabel('Number of epochs denoised at a time')
    ylabel('Number of pcs removed from data')
end

if saveFigures
    printnice(gcf, [1 300], fullfile(figureDir),'SuppFigure_Epochdur_Subjmean_vs_NPCS');
end