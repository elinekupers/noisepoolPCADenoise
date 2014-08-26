% define paths and data sets 
inputDataDir = '/Volumes/server/Projects/MEG/GLMdenoised/tmpmeg';
% noisepool selection by SNR, highpass filtered, 10 pcs removed
% bootstrapped 1000 x
fitDataStr   = 'b2fr_hpf2_fit10p1k_varyEpochs';
whichfun     = 1; %1 = broadband, 2 = stimulus-locked, 3=broadbandlog
condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
sessionNums = [11:12,3:6,9:10];

figuredir   = 'manuscript_figs/figure_controls';
savefigures = false;

%% Plot difference in SNR (post-pre) as a function of denoising epoch duration - Fig.11 A

epochDurs = [1,3,6,12,24,36,72,1080];
snr_diff = zeros(length(sessionNums),length(epochDurs),3);
nepochs  = zeros(1,length(sessionNums));

for k = 1:length(sessionNums)
    sessionDir = megGetDataPaths(sessionNums(k));
    % load fit file. this may take a while 
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile);

    % get all the results
    results_all = catcell(1,allResults);
    
    % get top 10 channels and use the same 10 channels for all epoch
    % durations
    pcchan = getTop10(results_all(1),whichfun);
    
    % total number of epochs for this session 
    nepochs(k) = length(results_all(1).opt.epochGroup);
    
    % compute the difference between pre and post
    for nn = 1:length(epochDurs)
        % get top 10 channels. note the choice of top 10 changes for each glm
        % result returned (for each epoch duration)
        %pcchan = getTop10(results_all(nn),whichfun);
        for icond = 1:3
            snr_pre  = getsignalnoise(results_all(nn).origmodel(whichfun),icond);
            snr_post = getsignalnoise(results_all(nn).finalmodel(whichfun),icond);
            snr_diff(k,nn,icond) = mean(snr_post(pcchan)-snr_pre(pcchan));
        end
    end
end

%%
% define colors - vary saturation across individual subjects
satValues = linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);

% set up figure and plot 
fH = figure('position',[0,300,250,600]);
for icond = 1:3 % for each condition 
    subplot(3,1,icond); hold on;
    % plot snr difference versus epoch duration 
    for nn = 1:length(sessionNums)
        plot([epochDurs(1:end-1),nepochs(nn)], snr_diff(nn,:,icond),'o-','color',squeeze(colorRGB(icond,nn,:)),'linewidth',2);
    end
    % format axes and make figure pretty 
    set(gca,'xscale','log');
    xlim([0.5,1500]); set(gca,'xtick',epochDurs,'xscale','log');
    ylim([-2,7]);
    ylabel('Difference in SNR (post-pre)');
    makeprettyaxes(gca,9,9);
end

if savefigures
    figurewrite(fullfile(figuredir,'figure_epochdur'),[],0,'.',1);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the surface of SNR difference as a function of epoch duration and
%% number of PCs removed. Fig. 11B 
% noisepool selection by SNR, highpass filtered, 10 PCs removed
% bootstrapped 100 x. This varies along one additional dimension than
% above (number of PCs removed). Ended up bootstrapping only 100x because
% the files were getting too huge. 
% see hpc/denoisescripthpc_megVaryEpochDur.m

fitDataStr   = 'b2fr_hpf2_fitfull75p1k_varyEpochs';
whichfun     = 1; %1 = broadband, 2 = broadbandlog
epochDurs = [1,3,6,12,24,36,72,1080];
npcs      = [5,10:10:70];

snr_diff = zeros(length(sessionNums),length(epochDurs),length(npcs),3);

for k = 1:length(sessionNums)
    sessionDir = megGetDataPaths(sessionNums(k));
    % load fit file. this may take a while
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile);

    % get all the results for a particular number of pc's removed
    for jj = 1:length(npcs)
        results_all = catcell(1,allResults(:,jj));
        
        % get top 10 channels
        pcchan = getTop10(results_all(1),whichfun);
        
        % compute the difference between pre and post
        for nn = 1:length(epochDurs)
            %pcchan = getTop10(results_all(nn),whichfun);
            for icond = 1:3               
                snr_pre  = getsignalnoise(results_all(nn).origmodel(whichfun),icond);
                snr_post = getsignalnoise(results_all(nn).finalmodel(whichfun),icond);
                snr_diff(k,nn,jj,icond) = mean(snr_post(pcchan)-snr_pre(pcchan));
            end
        end
    end
end

%% Plot 
fH = figure('position',[0,300,300,600]);
epochDurs = [1,3,6,12,24,36,72,1080];
npcs      = [5,10:10:70];
clims = [[0,4];[0,2];[0,2]];
% plot each condition as a separate panel 
for icond = 1:3
    subplot(3,1,icond);
    % plot the mean across subject as one square 
    sessionAvg = squeeze(mean(snr_diff(:,:,:,icond),1));
    imagesc(sessionAvg',clims(icond,:)); 
    % set up axes and format figure
    set(gca,'ydir','normal');
    makeprettyaxes(gca,9,9);
    set(gca,'xtick',1:length(epochDurs),'ytick',1:length(npcs),...
        'xticklabel',cellstr(num2str(epochDurs','%d')),'yticklabel',cellstr(num2str(npcs','%d')));
    %axis image; 
    ch = colorbar; makeprettyaxes(ch,9,9);
end

if savefigures
    figurewrite(fullfile(figuredir,'figure_epochdur_subjmean'),[],0,'.',1);
end