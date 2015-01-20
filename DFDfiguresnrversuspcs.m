function DFDfiguresnrversuspcs(sessionNums, conditionNumbers, plotBb, inputDataDir, whichFun, figureDir, saveFigures)

%% define paths and data sets

if notDefined('opt'),                   opt = struct(); end
if ~isfield(opt,'sessionNums'),         opt.sessionNums      = 1:8; end % Do all subjects
if ~isfield(opt,'conditionNumbers'),    opt.conditionNumbers = 1:6; end % Do all conditions
if ~isfield(opt,'plotBb'),              opt.plotBb           = true; end % true for broadband data, false for stimlocked data
if ~isfield(opt,'inputDataDir'),        opt.inputDataDir     = fullfile(DFDrootpath, 'data'); end % 
if ~isfield(opt,'whichFun'),            opt.whichFun         = 1; end % 
% if ~isfield(opt,'fitDataStr'),          opt.fitDataStr      = 'b2fr_hpf2_fit10p1k'; end % We need to automate this / link it to the opt. 
if ~isfield(opt,'figureDir'),           opt.figureDir        = fullfile(DFDrootpath,'figures'); end
if ~isfield(opt,'saveFigures'),         opt.saveFigures      = false; end

% fit data file string
if plotBb
    % noisepool selection by SNR, highpass filtered, 10 pcs removed
    % bootstrapped 1000 x. Broadband as evalfun
    fitDataStr   = 'b2fr_hpf2_fitfull75p1k';
%     fitDataStr   = 'b2fr_hpf2_fit10p1k';  % <------ FIX THIS HARDCODED PART
else
    % Same as above, but with stimulus locked as evalfun, and not hpf'ed
    fitDataStr   = 'b2frSL_fitfull75p1k'; % <------ FIX THIS HARDCODED PART
end

condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
axmax = 10; % how far to go out on the x-axis


%% SNR increase as a function of number of PCs removed, 3 example sessions - Fig. 6A
exampleSessions = [3,4,5];
linecolors = copper(157);

for k = 1:length(exampleSessions)
    % get session
    sessionDir = DFDgetdatapaths(exampleSessions(k),conditionNumbers,inputDataDir);
    % load fit file
    thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results','evalout');
    
    % get snr
    snr = abs(cat(3,evalout(:,whichFun).beta_md)) ./ cat(3,evalout(:,whichFun).beta_se);
    
    % plot for each condition
    for icond = 1:3
        subplot(4,3,(k-1)*3+icond); hold on;
        this_snr = squeeze(snr(icond,:,:))';
        % plot each channel's snr as a function of number of pc's
        for ic = 1:size(this_snr,2)
            plot(0:axmax,this_snr(1:axmax+1,ic),'color',linecolors(ic,:));
        end
        % plot snr change for top10 channels
        xvaltrend = mean(this_snr(:,results.pcchan{whichFun}),2);
        plot(0:axmax, xvaltrend(1:axmax+1,:), 'color', condColors(icond,:), 'linewidth',2);
        %plot(axmax+1, xvaltrend(51,:), 'o', 'color', condColors(icond,:));
        axis square; xlim([0,axmax]);
        if plotBb, ylim([0,15]); else ylim([0,50]); end
        makeprettyaxes(gca,9,9);
    end
end
if saveFigures
    if plotBb
        figurewrite(fullfile(figureDir,'Figure6SNRvPCsExampleSubjects'),[],0,'.',1);
    else 
        figurewrite(fullfile(figureDir,'Figure11SNRvPCsExampleSubjects_SL'),[],0,'.',1);
    end
end
%% SNR increase as a function of number of PCs removed for all subjects -
%% Fig. 6B

% get the trend for the top 10 channels of all sessions
% files might take a while to load!
snr_top10 = [];
for k = 1:length(sessionNums)
    sessionDir = DFDgetdatapaths(sessionNums(k),conditionNumbers,inputDataDir);
    % load fit file
    thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results','evalout');
    
    snr = abs(cat(3,evalout(:,whichFun).beta_md)) ./ cat(3,evalout(:,whichFun).beta_se);
    
    xvaltrend = [];
    for icond = 1:3
        this_snr = squeeze(snr(icond,:,1:11))';
        xvaltrend = cat(2, xvaltrend, mean(this_snr(:,results.pcchan{whichFun}),2));
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
    subplot(4,3,9+icond);hold on;
    for nn = 1:8 % for each subject
        plot(0:axmax, squeeze(snr_top10(:,icond,nn)), 'color', squeeze(colorRGB(icond,nn,:)));
    end
    axis square; xlim([0,axmax]);
    if plotBb
        ylim([0,12]); set(gca,'ytick',0:5:10);
    else
        ylim([0,40]); set(gca,'ytick',0:10:40);
    end
    title(ttls{icond});
    makeprettyaxes(gca,9,9);
end

if saveFigures
    if plotBb
        figurewrite(fullfile(figureDir,'Figure6BAllSubjs_sat_BB'),[],0,'.',1);
    else
        figurewrite(fullfile(figureDir,'Figure11BAllSubjs_sat_stimlocked'),[],0,'.',1);
    end
end

return