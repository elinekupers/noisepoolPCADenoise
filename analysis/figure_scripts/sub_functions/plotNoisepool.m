function fH = plotNoisepool(whichSubjects)

%% Function to plot the spatial map of channels that are in the noise pool

% fH = plotNoisePool(whichSubjects)

% AUTHORS. TITLE. JOURNAL. YEAR.

% This figure will show the channels that are in the noisepool in white
% (value of 1) and channels that are not in the noisepool as black (value
% of 0).

% INPUTS:
% whichSubjects   : Subjects you want to plot

% OUTPUTS:
% fH              : Figure handle

% Example:
% plotNoisePool(1:8)

%% Choices to make:                                              
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?

fH = figure('position',[1,600,1400,800]); set(gcf, 'Color','w');

%% Load denoised data of example subject
numrows = ceil(sqrt(length(whichSubjects)));
numcols = ceil(length(whichSubjects) / numrows);

for whichSubject = whichSubjects
    [data] = prepareData(dataDir,whichSubject,5);
    bb = data{1};

    % Get noisepool data
    noisepool = double(bb.results.noisepool);
    noisepool = to157chan(noisepool,~bb.badChannels,'zeros');

    %% Plot it
    % Define binary color range
    clims = [0,1];
    
    subplot(numrows, numcols, find(whichSubject==whichSubjects));
    % Plot
    [~,ch] = megPlotMap(noisepool,clims,gcf,'gray', ...
        sprintf('S%d',whichSubject), [], [], 'interpmethod', 'nearest');
    makeprettyaxes(ch,9,9);
    colorbar('YTick',[0 1], 'YTicklabel',{'Not in noise pool','In noise pool'}, 'FontSize', 9)

    %% Save it
    figurewrite(sprintf(fullfile(figureDir,'noisepool_subject%d'),whichSubject),[],0,'.',1);

end