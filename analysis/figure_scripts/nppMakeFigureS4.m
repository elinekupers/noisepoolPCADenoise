function nppMakeFigureS4()
%% Function to reproduce supplementary Figure 4 SNR pre-post denoising
% when removing pcs from a synthetic data set that contains 10 bases
% functions
%
% nppMakeFigureS4()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will show the simulated broadband response before and after denoising.
% Bootstrapped fullfield signal (mean across bootstraps) and noise component (std across bootstraps)
% and the difference between the two distributions are plotted before and
% after denoising.
%
% This function assumes that data is downloaded with the nppdownloaddata
% function.

%% Choices to make:
whichSubject    = 99;     % Subject 99 has the synthetic dataset.
figureDir       = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
ylimsSF1D       = [-1, 8];

% Define plotting parameters
colors          = nppGetColors(4);
axmax           = 75;    % How far out do you want to plot the number of PCs

fprintf('Load data subject %d \n', whichSubject);
% Load data, design, and get example subject
[data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,4);

% Define conditions: Full, right, left, off
condEpochs1 = {design{1}(:,1)==1, design{1}(:,2)==1, design{1}(:,3)==1, all(design{1}==0,2)};
condEpochs2 = {design{2}(:,1)==1, design{2}(:,2)==1, design{2}(:,3)==1, all(design{2}==0,2)};

condEpochs = {condEpochs1 condEpochs2};

%% Create power spectrum, before denoising

% Set up figure
fH(1) = figure('position',[0,300,500,500]); clf(fH(1)); set(fH(1), 'Name', 'S4A, Power spectrum', 'NumberTitle', 'off');

% Define axes
f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];
yt = -4:2;
yl=[yt(1),yt(end)];

% Compute spectrum
spec = abs(fft(squeeze(data{1}(exampleIndex,:,:))))/size(data{1},2)*2;

hold on;
for ii = [1,4] % For both-field and blank stimulus
    
    % Compute power for epochs corresponding to a condition and
    % trim data to specific frequencies
    this_data = spec(:,condEpochs{1}{ii}).^2;
    this_data = this_data(fok+1,:);
    
    % Compute median and confidence interval across epochs
    mn = prctile(this_data,[16,50,84],2);
    
    % Plot median
    plot(fok, mn(:,2),  '-',  'Color', colors(ii,:), 'LineWidth', 1);
    
end

% Format x and y axes
set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
set(gca, 'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');

% Label figure, add stimulus harmonic lines, and make it look nice
xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
title('Before denoising');
yl2 = get(gca, 'YLim');
for ii =12:12:180, plot([ii ii], yl2, 'k--'); end
makeprettyaxes(gca,9,9);

if saveFigures
    figurewrite(sprintf(fullfile(figureDir,'figureS4AFullSpectrumChannel%d'),exampleChannel),[],0,'.',1);
end

%% Broadband before and after denoising

% Set up figure
fH(2) = figure; set(fH,'position',[0,300,200,350], 'Name', 'S4B and C, spectra zoom before and after denoising', 'NumberTitle', 'off')

% Rescale axes and remove harmonics
xl = [60 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 | mod(f,72) < 2 | mod(f,96) < 2 | mod(f,108) < 2 | mod(f,144)<2) = [];
xt = [];
yt = -5:-1;
yl=[yt(1),yt(end)];

for dd = 1:2
    subplot(2,1,dd);
    
    % compute spectrum
    spec = abs(fft(squeeze(data{dd}(exampleIndex,:,:))))/size(data{dd},2)*2;
    
    hold on;
    for conditions = [1,4]
        % compute power for epochs corresponding to a condition and
        % trim data to specific frequencies
        this_data = spec(:,condEpochs{dd}{conditions}).^2;
        this_data = this_data(fok+1,:);
        
        mn = prctile(this_data,[16,50,84],2);
        
        % Plot the data
        plot(fok, mn(:,2),  '-',  'Color', colors(conditions,:), 'LineWidth', 1);
        set(gca,'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');
        set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
        
        % label figure, add stimulus harmonic lines, and make it look nice
        xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
        if dd == 1, title('Before denoising'), else title('After denoising'); end
        yl2 = get(gca, 'YLim');
        for ii =12:12:180, plot([ii ii], yl2, 'k--'); end
        makeprettyaxes(gca,9,9);
        
    end
end

if saveFigures
    figurewrite(sprintf(fullfile(figureDir,'figureS1BBeforeAfterSpectrumChannel%d'),exampleChannel),[],0,'.',1);
end

%% Plot SNR vs number of PCs change for all channels
data = prepareData(dataDir,whichSubject,7);
fH(3) = plotSNRvsPCsAllSubjectsPanel7B({data},colors,axmax,figureDir,saveFigures,ylimsSF1D, 'Figure SF4D');

%% Plot noise pool
fH(4) = plotNoisepool(whichSubject);

%% Plot topographic maps before and after denoising
fH(5) = nppMakeFigure9(whichSubject); 

