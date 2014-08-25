%% define paths and data sets
plotbb = true; % true for broadband data, false for stimlocked data
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
% fit data file string
if plotbb
    fitDataStr   = 'b2fr_hpf2_fitfull75p1k';
else
    fitDataStr   = 'b2frSL_fitfull75p1k';
end
whichfun     = 1;
condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
sessionNums  = [11,12, 3:6, 9:10];%[1:6,9,10];
axmax = 10;

figuredir = 'manuscript_figs/figure_snrVpcs';

%% SNR increase as a function of number of PCs removed, 3 example sessions - Fig. 6A
exampleSessions = [5,6,9];
linecolors = copper(157);

for k = 1:length(exampleSessions)
    % get session
    sessionDir = megGetDataPaths(exampleSessions(k));
    % load fit file
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results','evalout');
    
    % get snr
    snr = abs(cat(3,evalout(:,whichfun).beta_md)) ./ cat(3,evalout(:,whichfun).beta_se);
    
    % plot for each condition
    for icond = 1:3
        subplot(3,3,(k-1)*3+icond); hold on;
        this_snr = squeeze(snr(icond,:,:))';
        % plot each channel's snr as a function of number of pc's
        for ic = 1:size(this_snr,2)
            plot(0:axmax,this_snr(1:axmax+1,ic),'color',linecolors(ic,:));
        end
        % plot snr change for top10 channels
        xvaltrend = mean(this_snr(:,results.pcchan{whichfun}),2);
        plot(0:axmax, xvaltrend(1:axmax+1,:), 'color', condColors(icond,:), 'linewidth',2);
        %plot(axmax+1, xvaltrend(51,:), 'o', 'color', condColors(icond,:));
        axis square; xlim([0,axmax]);
        if plotbb, ylim([0,15]); else ylim([0,50]); end
        makeprettyaxes(gca,9,9);
    end
end

%figurewrite(fullfile(figuredir,'figure_examples'),[],0,'.',1);

%% SNR increase as a function of number of PCs removed for all subjects -
%% Fig. 6B

% get the trend for the top 10 channels of all sessions
% files might take a while to load!
snr_top10 = [];
for k = 1:length(sessionNums)
    sessionDir = megGetDataPaths(sessionNums(k));
    % load fit file
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results','evalout');
    
    snr = abs(cat(3,evalout(:,whichfun).beta_md)) ./ cat(3,evalout(:,whichfun).beta_se);
    
    xvaltrend = [];
    for icond = 1:3
        this_snr = squeeze(snr(icond,:,:))';
        xvaltrend = cat(2, xvaltrend, mean(this_snr(:,results.pcchan{whichfun}),2));
    end
    snr_top10 = cat(3,snr_top10,xvaltrend);
end

%% Plot them together

% define colors - vary saturation for different subjects
satValues = linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);
ttls = {'FULL','RIGHT','LEFT'};
% plot for each condition
for icond = 1:3
    subplot(1,3,icond);hold on;
    for nn = 1:8 % for each subject
        plot(0:75, squeeze(snr_top10(:,icond,nn)), 'color', squeeze(colorRGB(icond,nn,:)));
    end
    axis square; xlim([0,axmax]);
    if plotbb
        ylim([0,12]); set(gca,'ytick',0:5:10);
    else
        ylim([0,40]); set(gca,'ytick',0:10:40);
    end
    title(ttls{icond});
    makeprettyaxes(gca,9,9);
end
%figurewrite(fullfile(figuredir,'figure_allsubjs_sat'),[],0,'.',1);
