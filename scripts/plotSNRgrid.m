function plotSNRgrid(allvals,npools,npcs, ...
    whichfun, doTop10, plotType, whichConds, funXchan)

% plot stuff
figure('position',[200,200,600,600]); 
imagesc(1:length(npcs),1:length(npools),allvals); 
xlabel('Number of PCs removed'); 
ylabel('Number of Channels in Noise pool');
if length(whichConds) > 1, condName = 'all'; else condName = num2str(whichConds); end
if doTop10, st = 'Top 10'; else st = 'Non Noise'; end
ttl = {sprintf('Diff in %s (post - pre)',upper(plotType));...
    sprintf('Cond: %s, %s across %s channels', condName, func2str(funXchan), st)};
title(ttl);
set(gca,'ydir','normal');
makeprettyaxes(gca,14,14); 
set(gca,'xtick',1:length(npcs),'ytick',1:length(npools),...
    'xticklabel',cellstr(num2str(npcs','%d')),'yticklabel',cellstr(num2str(npools','%d')));
axis image; colorbar;