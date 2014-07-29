clear all;
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
outputFigDir = 'megfigs';
sessionNums  = 1:10;
fitDataStr    = 'b2fr_epochGroup6so_fitall';
whichfun      = 3;        % which fit (usually only 1)
doTop10       = true;
plotType      = 'n';
whichConds    = 1:3;
funXchan      = @mean;

evalfuncs = {'bb','sl','bblog'};
printFigsToFile = true;

%%
allvals2 = [];
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    sessionDir = megGetDataPaths(sessionNums(k));
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile);
    
    allvals = getSNRgrid(allresults,npools,npcs,...
        whichfun,doTop10,plotType,whichConds,funXchan);
    
    allvals2 = cat(3,allvals2,allvals);
    % write file
    if printFigsToFile
        savename = sprintf('%s_%s_%02d_%s%s', plotType, evalfuncs{whichfun}, sessionNums(k), sessionDir, fitDataStr);
        figurewrite(savename, [],[], outputFigDir, 0);
    else
        pause;
    end
end

%%
plotSNRgrid(mean(allvals2,3),npools,npcs,...
    whichfun,doTop10,plotType,whichConds,funXchan);

if printFigsToFile
    savename = sprintf('%s_%s_%s_meanAll', plotType, evalfuncs{whichfun}, fitDataStr);
    figurewrite(savename, [],[], outputFigDir, 0);
end