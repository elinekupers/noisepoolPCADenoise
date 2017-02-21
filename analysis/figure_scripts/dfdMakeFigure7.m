function dfdMakeFigure7()
%% Function to reproduce Figure 7AB S, N pre-post denoising
% for top ten channels of all subjects
%
% dfdMakeFigure7()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% This figure will show subject's the broadband spectra before and after denoising from 60-150 Hz.
% Bootstrapped fullfield signal (mean across bootstraps) and noise component (std across bootstraps) 
% and the difference between the two distributions are plotted before and
% after denoising. 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:
whichSubjects    = [1:8];     % Subject 1 has the example channel.
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figureNumber    = 7;
                                         
% Define plotting parameters
colors          = dfdGetColors(3);

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
fH = figure('position',[0,300,500,400]); set(gcf, 'Color','w', 'Name', 'Figure 7A&B', 'NumberTitle', 'off');
datatypes = {'Noise','Signal'};
for t = 1:numel(datatypes);
    for icond = 1:3
        subplot(numel(datatypes),3,((t-1)*3+icond))
        plotBeforeAfter(allresults,1,allpcchan,datatypes{t},icond,[],squeeze(colorRGB(icond,:,:)));
        xlim([0.5,2.5]);
        makeprettyaxes(gca,9,9);
        ylim([-2,10]);
    end
end

if saveFigures
        figurewrite(fullfile(figureDir,'Figure7_s_n_full_sat'),[],0,'.',1);
end


