function fH = plotEpochLengthVersusSNR(whichSubjects,dataAll,epochDurs,condColors,saveFigures,figureDir)


for npsIdx = 2; % Corresponding to 10 pcs (1 = 5 PCs, 2 = 10 PCs, in steps of 10 up to 70)
    snr_diff = NaN(length(whichSubjects),length(epochDurs),3);
    nepochs  = NaN(1,length(whichSubjects));
    
    for whichSubject = whichSubjects
        
        % get all the results
        results_all = catcell(1,dataAll{whichSubject}{1}.allResults(:,npsIdx));
        npcs = dataAll{1}{1}.npcs;
        
        % get top 10 channels and use the same 10 channels for all epoch
        % durations
        pcchan = getTop10(results_all(1));
        
        % total number of epochs for this session
        nepochs(whichSubject) = length(results_all(1).opt.epochgroup);
        
        % compute the difference between pre and post
        for nn = 1:length(epochDurs)
            
            bb_signal_pre = results_all(nn).origmodel(1).beta_md(:,pcchan);
            bb_noise_pre  = results_all(nn).origmodel(1).beta_se(:,pcchan);
            
            bb_signal_post = results_all(nn).finalmodel(1).beta_md(:,pcchan);
            bb_noise_post  = results_all(nn).finalmodel(1).beta_se(:,pcchan);
            
            bb_snr_pre     = bb_signal_pre./bb_noise_pre;
            bb_snr_post    = bb_signal_post./bb_noise_post;
             
            snr_diff(whichSubject,nn,:) = mean((bb_snr_post-bb_snr_pre),2);
            
        end
        
    end
    
    clear allResults
    %%
    % define colors - vary saturation across individual subjects
    satValues = 1-linspace(0.1,1,8);
    colorRGB = varysat(condColors,satValues);
    ylims = {[0,7],[0,5],[0,5]};
    % set up figure and plot
    fH = figure('position',[0,300,250,600]); clf;
    set(fH, 'Name', sprintf('%d PCs used for denoising', npcs(npsIdx)), 'NumberTitle', 'off');
    for icond = 1:3 % for each condition
        subplot(3,1,icond); hold on;
        % plot snr difference versus epoch duration
        for nn = 1:length(whichSubjects)
            plot([epochDurs(1:end-1),nepochs(nn)], snr_diff(nn,:,icond),'o-','color',squeeze(colorRGB(icond,nn,:)),'linewidth',2);
        end
        
        plot([epochDurs(1:end-1),nepochs(nn)], mean(snr_diff(:,:,icond),1),'k-','linewidth',2);
        
        % format axes and make figure pretty
        set(gca,'xscale','log');
        xlim([0.5,1500]); set(gca,'xtick',epochDurs,'xscale','log');
        ylim(ylims{icond}); set(gca, 'YTick', [-2:2:ylims{icond}(2)])
        ylabel('Difference in SNR (post-pre)');
        xlabel('Number of epochs denoised at a time');
        axis square;
        makeprettyaxes(gca,9,9);
    end
    
    if saveFigures
        figurewrite(fullfile(figureDir,sprintf('figure_epochdur_npcsUsedToDenoise%d',npsIdx)),[],0,'.',1);
    end
end
