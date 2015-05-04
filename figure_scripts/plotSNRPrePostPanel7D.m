function fH = plotSNRPrePostPanel7D(dataAll,whichSubjects,condColors,figureDir,saveFigures)



% get results for everybody and top10 channels for everybody
for k = 1:length(whichSubjects)
    allpcchan{k} = getTop10(dataAll{k}{1}.results);
    allresults{k} = dataAll{k}{1}.results;
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
    ylim([0,12]); % for SL plot, take ylim([0,40])
end

if saveFigures
   figurewrite(fullfile(figureDir,'Figure7D_snrfull_sat'),[],0,'.',1);
end