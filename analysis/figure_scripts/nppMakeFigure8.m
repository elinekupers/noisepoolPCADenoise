function nppMakeFigure8()
%% Function to reproduce Figure 8AB S, N pre-post denoising
% for top ten channels of all subjects
%
% nppMakeFigure8()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will show subject's the broadband spectra before and after denoising from 60-150 Hz.
% Bootstrapped fullfield signal (mean across bootstraps) and noise component (std across bootstraps) 
% and the difference between the two distributions are plotted before and
% after denoising. 
%
% This function assumes that data is downloaded with the nppDownloaddata
% function. 

%% Choices to make:
whichSubjects   = 1:8;     % Subject 1 has the example channel.
figureDir       = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figureNumber    = 8;
if whichSubjects == 99; yl = [0,0.01]; else yl = [-2,10]; end
                                         
% Define plotting parameters
colors          = nppGetColors(3);

%% Get data
dataAll = [];
for whichSubject = whichSubjects
    fprintf('Load data subject %d \n', whichSubject);
    % Load data, design, and get example subject
    dataAll{whichSubject} = prepareData(dataDir,whichSubject,figureNumber);
end

%% Plot SNR vs number of PCs change for all channels 

% Get results for everybody and top10 channels for everybody
for k = whichSubjects
    allpcchan{k} = getTop10(dataAll{k}.results);
    allresults{k} = dataAll{k}.results;
end

% get colors for plotting
% vary saturation for different subjects
satValues = 1-linspace(0.1,1,8);
colorRGB = varysat(colors,satValues);

% plot before and after
fH = figure('position',[0,300,500,400]); set(gcf, 'Color','w', 'Name', 'Figure 8A&B', 'NumberTitle', 'off');
datatypes = {'Noise','Signal'};
for t = 1:numel(datatypes);
    for icond = 1:3
        subplot(numel(datatypes),3,((t-1)*3+icond))
        plotBeforeAfter(allresults,1,allpcchan,datatypes{t},icond,[],squeeze(colorRGB(icond,:,:)));
        xlim([0.5,2.5]);
        makeprettyaxes(gca,9,9);
        ylim(yl);
    end
end

if saveFigures
        figurewrite(fullfile(figureDir,'Figure8_s_n_full_sat'),[],0,'.',1);
end


