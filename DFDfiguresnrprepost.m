function DFDfiguresnrprepost(sessionNums, conditionNumbers, plotBb, inputDataDir, whichFun, figureDir, saveFigures)

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


if plotBb
    % noisepool selection by SNR, highpass filtered, 10 pcs removed
    % bootstrapped 1000 x. Broadband as evalfun
    fitDataStr = 'b2fr_hpf2_fit10p1k'; 
else
    % Same as above, but with stimulus locked as evalfun, and not highpass
    % filtered
    fitDataStr = 'b2frSL_fit10p1k';    
end

condColors   = [63, 121, 204; 228, 65, 69; 116,183,74]/255;

%% S, N, and SNR shown separately, before versus after denoising with 10
%% PCs. For 3 example sessions - Fig. 7A,B,C

exampleSessions = [3,4,5]; % Helena's [5,6,9]
for k = 1:length(exampleSessions)
    fprintf(' session %d \n', exampleSessions(k));
    sessionDir = DFDgetdatapaths(exampleSessions(k),conditionNumbers,inputDataDir);
    
    % load fit file
    thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results');
    
    % top 10 channels
    pcchan = getTop10(results,whichFun);
    %pcchan = results.pcchan{whichFun};
    %pcchan = ~results.noisepool;
    %pcchanfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,'b2fr_hpf2_fitfull75'));
    %tmp = load(pcchanfile); pcchan = tmp.results.pcchan{whichFun};
    
    % signal and noise before denoising
    ab_signal1 = abs(results.origmodel(whichFun).beta_md(:,pcchan));
    ab_noise1  = results.origmodel(whichFun).beta_se(:,pcchan);
    % signal and noise after denoising 
    ab_signal2 = abs(results.finalmodel(whichFun).beta_md(:,pcchan));
    ab_noise2  = results.finalmodel(whichFun).beta_se(:,pcchan);
    % snr = signal/denoise
    ab_snr1    = ab_signal1./ab_noise1;
    ab_snr2    = ab_signal2./ab_noise2;
    
    % plot each condition as a different color 
    % signal 
    subplot(3,length(exampleSessions),k); cla; hold on;
    for nn = 1:3
        plot(ab_signal1(nn,:),ab_signal2(nn,:),'o','color',condColors(nn,:));
    end
    axis square;
    axismax = max([ab_signal1(:);ab_signal2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    title(sprintf('S%d : signal', exampleSessions(k)));
    makeprettyaxes(gca,9,9);
    
    % noise
    subplot(3,length(exampleSessions),k+length(exampleSessions)); cla; hold on;
    for nn = 1:3
        plot(ab_noise1(nn,:),ab_noise2(nn,:),'o','color',condColors(nn,:));
    end
    axismax = max([ab_noise1(:); ab_noise2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : noise', exampleSessions(k)));
    makeprettyaxes(gca,9,9);
    
    % snr
    subplot(3,length(exampleSessions),k+2*length(exampleSessions)); cla; hold on;
    for nn = 1:3
        plot(ab_snr1(nn,:),ab_snr2(nn,:),'o','color',condColors(nn,:));
    end
    axismax = max([ab_snr1(:); ab_snr2(:)])*1.2;
    xlim([0,axismax]); ylim([0,axismax]); line([0,axismax],[0,axismax],'color','k');
    axis square;
    title(sprintf('S%d : SNR', exampleSessions(k)));
    makeprettyaxes(gca,9,9);
    
    drawnow;
end

if saveFigures
    if plotBb
        figurewrite(fullfile(figureDir,'Figure7abc_snrexamples'),[],0,'.',1);
    else
        figurewrite(fullfile(figureDir,'Figure12abc_snrexamples_SL'),[],0,'.',1);
    end
end

%% Plot changes in SNR before and after denoising, showing all sessions
%% together - Fig. 7D

% get results for everybody and top10 channels for everybody
clear allpcchan
for k = 1:length(sessionNums)
    sessionDir = DFDgetdatapaths(sessionNums(k),conditionNumbers,inputDataDir);
    % load fit file
    thisfile = fullfile(inputDataDir,'savedProcData',sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'results');
    allresults{k} = results;
    
    %pcchanfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,'b2fr_hpf2_fitfull75'));
    %tmp = load(pcchanfile); pcchan{k} = tmp.results.pcchan{whichFun};
    allpcchan{k} = getTop10(results,whichFun);
end

% get colors for plotting
% vary saturation for different subjects 
satValues = linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);

% plot before and after
fH = figure('position',[0,300,500,200]);
for icond = 1:3
    subplot(1,3,icond);
    plotBeforeAfter(allresults,1,allpcchan,'snr',icond,[],squeeze(colorRGB(icond,:,:)));
    xlim([0.5,2.5]);
    makeprettyaxes(gca,9,9);
    if plotBb, ylim([0,12]); else ylim([0,40]); end
end

if saveFigures
    if plotBb
        figurewrite(fullfile(figureDir,'Figure7D_snrfull_sat'),[],0,'.',1);
    else
        figurewrite(fullfile(figureDir,'Figure12D_snrfull_sat_SL'),[],0,'.',1);
    end
end

return
