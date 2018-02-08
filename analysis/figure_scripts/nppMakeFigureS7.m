function nppMakeFigureS7()
%% Function to reproduce Supplementary Figure 7AB S, N pre-post denoising
% for top ten channels of all subjects for stimulus locked signal
%
% nppMakeFigureS7()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will SNR, signal (mean across bootstraps) and noise (std
% across bootstraps) components of subject's stimulus.

%
% This function assumes that data is downloaded with the nppdownloaddata
% function. 

%% Choices to make:
whichSubjects    = [1:8];     % Subject 1 has the example channel.
figureDir       = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;  % Save figures in the figure folder?
figureNumber    = 'SF7';
                                         
% Define plotting parameters
colors          = nppGetColors(3);

%% Get data
dataAll = [];
for whichSubject = whichSubjects
    fprintf('Load data subject %d \n', whichSubject);
    % Load data, design, and get example subject
    dataAll = prepareData(dataDir,whichSubject,figureNumber); 

    % Get results for everybody and top10 channels for everybody
    allpcchan{whichSubject} = getTop10(dataAll.results);
    allresults{whichSubject} = dataAll.results;
end

clear dataAll

% get colors for plotting
% vary saturation for different subjects
satValues = 1-linspace(0.1,1,8);
colorRGB = varysat(colors,satValues);

% plot before and after
fH = figure('position',[0,300,400,600],'Name', 'Figure S7', 'NumberTitle', 'off'); set(gcf, 'Color','w');
datatypes = {'SNR','Signal','Noise'};
for t = 1:numel(datatypes);
    for icond = 1:3
        subplot(numel(datatypes),3,((t-1)*3+icond))
        plotBeforeAfter(allresults,1,allpcchan,datatypes{t},icond,[],squeeze(colorRGB(icond,:,:)));
        xlim([0.5,2.5]);
        makeprettyaxes(gca,9,9);
        if t==1; yt = [0,40]; elseif t==2; yt= [0,130]; else yt = [0,7]; end
        ylim(yt);
    end
end

if saveFigures
   figurewrite(fullfile(figureDir,'FigureS7_s_n_full_sat'),[],0,'.',1);
end


