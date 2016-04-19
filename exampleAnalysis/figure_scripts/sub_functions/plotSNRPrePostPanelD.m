function fH = plotSNRPrePostPanelD(dataAll,whichSubjects,condColors,figureDir,saveFigures,figureNumber)

% get results for everybody and top10 channels for everybody
for k = whichSubjects
    allpcchan{k} = getTop10(dataAll{k}{1}.results);
    allresults{k} = dataAll{k}{1}.results;
end

% get colors for plotting
% vary saturation for different subjects
satValues = 1-linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);

% plot before and after
fH = figure('position',[0,300,500,200]);
    for icond = 1:3
        subplot(1,3,icond);
        plotBeforeAfter(allresults,1,allpcchan,'SNR',icond,[],squeeze(colorRGB(icond,:,:)));
        xlim([0.5,2.5]);
        makeprettyaxes(gca,9,9);
        if isequal(figureNumber,7);  ylim([-2,10]);
        else ylim([0,40]); end;
    end

if saveFigures
    if isequal(figureNumber,7)
        figurewrite(fullfile(figureDir,'Figure7D_snr_full_sat'),[],0,'.',1);
    elseif isequal(figureNumber,11)
        figurewrite(fullfile(figureDir,'Figure11D_snr_full_sat'),[],0,'.',1);
    end
end