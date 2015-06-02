function dfdMakeFigure12()
%% Function to reproduce Figure 12A SNR against PCs removed for 3 example subjects 
% and Figure 12b, SNR for top ten channels against PCs removed for all
% subjects in stimulus-locked signals
%
% dfdMakeFigure12()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show 2 ways of plotting SNR values against the number of PCs
% removed from the data. In Figure 12A, SNR values of all channels from the 
% three example subjects (corresponding to 3,4,5) will be plotted against
% the number of PCs removed. The median is plotted in the color of the
% condition. This function needs the full dataset
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:
whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'figures');% Where to save images?
saveFigures          = true;     % Save figures in the figure folder?
exampleSessions      = [3,4,5];  % Helena's plot contained subjects [5,6,9]
condColors           = [63, 121, 204; 228, 65, 69; 116,183,74]/255;
linecolors = copper(157);
dataAll              = [];

%% Load data all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,11);  
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
     
end

%% SNR increase as a function of number of PCs for 3 example subjects
for k = exampleSessions
      
    data = dataAll{k};
    results = data{1}.results;
    
    % get snr
    snr = abs(cat(3,results.evalout.beta_md)) ./ cat(3,results.evalout.beta_se);
    
    % plot for each condition
    for icond = 1:3
        subplot(8,3,(k-1)*3+icond); hold on;
        this_snr = squeeze(snr(icond,:,:))';
        % plot each channel's snr as a function of number of pc's
        for ic = 1:size(this_snr,2)
            plot(0:axmax,this_snr(1:axmax+1,ic),'color',linecolors(ic,:));
        end
        % plot snr change for top10 channels
        xvaltrend = mean(this_snr(:,results.pcchan),2);
        plot(0:axmax, xvaltrend(1:axmax+1,:), 'color', condColors(icond,:), 'linewidth',2);
        %plot(axmax+1, xvaltrend(51,:), 'o', 'color', condColors(icond,:));
        axis square; xlim([0,axmax]);
        ylim([0,50]);
        makeprettyaxes(gca,9,9);
    end
end

if saveFigures
   figurewrite(fullfile(figureDir,'Figure12SNRvPCsExampleSubjects_SL'),[],0,'.',1);
end

%% SNR increase as a function of number of PCs removed for all subjects -
%% Fig. 12B

% get the trend for the top 10 channels of all sessions
% files might take a while to load!
snr_top10 = [];
for k = 1:length(whichSubjects)
    
    data = dataAll{k};
    results = data{1}.results;
    
    snr = abs(cat(3,results.evalout.beta_md)) ./ cat(3,results.evalout.beta_se);
    
    xvaltrend = [];
    for icond = 1:3
        this_snr = squeeze(snr(icond,:,1:11))';
        xvaltrend = cat(2, xvaltrend, mean(this_snr(:,results.pcchan),2));
    end
    snr_top10 = cat(3,snr_top10,xvaltrend);
end

% define colors - vary saturation for different subjects
satValues = linspace(0.1,1,8);
colorRGB = varysat(condColors,satValues);
ttls = {'FULL','RIGHT','LEFT'};

% plot for each condition
for icond = 1:3
    subplot(1,3,icond);hold on;
    for nn = 1:8 % for each subject
        plot(0:axmax, squeeze(snr_top10(:,icond,nn)), 'color', squeeze(colorRGB(icond,nn,:)));
    end
    axis square; xlim([0,axmax]);
    ylim([0,40]); set(gca,'ytick',0:10:40);
    title(ttls{icond});
    makeprettyaxes(gca,9,9);
end

if saveFigures
    figurewrite(fullfile(figureDir,'Figure12BSNRvPCsAllSubjects_SL'),[],0,'.',1);
end

return
