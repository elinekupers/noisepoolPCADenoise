% Plots the improvement in SNR as a function of the number of channels in
% the noise pool and the number of PCs removed 

clear all;
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
outputFigDir = 'megfigs';
sessionNums  = 11:12;%[1:6,9:10];
fitDataStr    = 'b2fr_hpf2_fitall';
whichfun      = 1;        % which evafun (usually only 1)
doTop10       = true;
plotType      = 'snr';
whichConds    = 1:3;
funXchan      = @mean;

evalfuncs = {'bb','sl','bblog'};
printFigsToFile = true;

%% Plot individual sessions 
allvals2 = [];
for k = 1:length(sessionNums)
    fprintf(' session %d \n', sessionNums(k));
    sessionDir = megGetDataPaths(sessionNums(k));
    % load data 
    thisfile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    disp(thisfile); load(thisfile,'allresults','npools','npcs');
    
    % put the results on a grid and plot it 
    allvals = getSNRgrid(allresults,npools,npcs,...
        whichfun,doTop10,plotType,whichConds,funXchan);
    
    % concatenate across sessions
    allvals2 = cat(3,allvals2,allvals);
    
    % write out figure 
    if printFigsToFile
        savename = sprintf('%s_%s_%02d_%s%s', plotType, evalfuncs{whichfun}, sessionNums(k), sessionDir, fitDataStr);
        figurewrite(savename, [],[], outputFigDir, 0);
    else
        pause;
    end
end

%% Plot the mean across sessions
plotSNRgrid(mean(allvals2,3),npools,npcs,...
    whichfun,doTop10,plotType,whichConds,funXchan);

if printFigsToFile
    savename = sprintf('%s_%s_%s_meanAll', plotType, evalfuncs{whichfun}, fitDataStr);
    figurewrite(savename, [],[], outputFigDir, 0);
end