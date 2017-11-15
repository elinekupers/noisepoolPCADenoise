function dfdMakeFigureS1()
%% Function to reproduce supplementary Figure 3 SNR pre-post denoising
% when removing pcs from a synthetic data set that contains 5 or 10 bases
% functions
%
% dfdMakeFigureS1()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show the simulated broadband response before and after denoising.
% Bootstrapped fullfield signal (mean across bootstraps) and noise component (std across bootstraps)
% and the difference between the two distributions are plotted before and
% after denoising.
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function.

%% Choices to make:
whichSubject    = 99;     % Subject 99 has the synthetic dataset.
exampleSessions = 1;      % this session contains 10 basis functions
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figureNumber    = 'SF1';
ylimsSF1D       = [-1, 5]; % for pink noise data we have to lower the ylimita (was [-20 20]);

% Define plotting parameters
colors          = dfdGetColors(4);
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
fH(1) = figure('position',[0,300,500,500]); clf(fH(1)); set(fH(1), 'Name', 'S1A, Power spectrum', 'NumberTitle', 'off');

% Define axes
f = (0:999);
xl = [8 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 ) = [];
xt = [12:12:72, 96,144];
yt = -2:2;
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
    figurewrite(sprintf(fullfile(figureDir,'figureS1AFullSpectrumChannel%d'),exampleChannel),[],0,'.',1);
end

%% Broadband before and after denoising

% Get data
[data,~,exampleIndex] = prepareData(dataDir,whichSubject,6);

% Set up figure
fH(2) = figure; set(fH,'position',[0,300,200,350], 'Name', 'S1B and C, spectra zoom before and after denoising', 'NumberTitle', 'off')

% Rescale axes and remove harmonics
xl = [60 150];
fok = f;
fok(f<=xl(1) | f>=xl(2) | mod(f,60) < 2 | mod(f,60) > 58 | mod(f,72) < 2 | mod(f,96) < 2 | mod(f,108) < 2 | mod(f,144)<2) = [];
xt = [];
yt = -3:0;
yl=[yt(1),yt(end)];

for dd = 1:2
    subplot(2,1,dd);
    
    % compute spectrum
    spec = abs(fft(squeeze(data{dd}(exampleIndex,:,:))))/size(data{dd},2)*2;
    
    hold on;
    for ii = [1,4]
        % compute power for epochs corresponding to a condition and
        % trim data to specific frequencies
        this_data = spec(:,condEpochs{dd}{ii}).^2;
        this_data = this_data(fok+1,:);
        
        mn = prctile(this_data,[16,50,84],2);
        
        % Plot the data
        plot(fok, mn(:,2),  '-',  'Color', colors(ii,:), 'LineWidth', 1);
        set(gca,'ytick',10.^yt, 'ylim',10.^yl,'YScale', 'log');
        set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'log', 'YScale','log');
        
        % label figure, add stimulus harmonic lines, and make it look nice
        xlabel('Frequency (Hz)'); ylabel('Power (fT^2)');
        title('After denoising');
        yl2 = get(gca, 'YLim');
        for ii =12:12:180, plot([ii ii], yl2, 'k--'); end
        makeprettyaxes(gca,9,9);
        
    end
end

%% Plot SNR vs number of PCs change for all channels
data = prepareData(dataDir,whichSubject,'SF1');
fH(3) = plotSNRvsPCsAllSubjectsPanel7B({data},colors,axmax,figureDir,saveFigures,ylimsSF1D, 'Figure SF1D');

%% Plot noise pool
fH(4) = plotNoisepool(whichSubject);

%% Plot topographic maps before and after denoising
fH(5) = dfdMakeFigure9(whichSubject); for ax = 2:2:16; set(fH(5).Children(ax),'CLim', [-20 20]); end
                                        for ax = 1:2:15; set(fH(5).Children(ax), 'YTick', [-20 -10 0 10 20]); end

