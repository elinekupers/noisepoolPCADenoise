function dfdMakeFigure()
%% Function to reproduce Figure 9 (Spatialmap) for right minus left denoised activations 
%
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the broadband signals in the contrast right minus left
% after using the denoising algorithm. 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:

whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'data');
figureDir            = fullfile(dfdRootPath, 'figures');
saveFigures          = true;   % Save figures in the figure folder?
dataAll              = [];

%% Check whether we got our preprocessed data matrices

for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,9);
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
end                                                       

%% Make spatial map figures:
% Plot right minus left stimulation, SL and BB pre and post denosing separately)
figure('position',[1,600,1400,800]);
whichmodel = 'finalmodel';
for k = 1:length(whichSubjects)

    data = dataAll{k};
    results = data{1}.results;
    
    subplot(2,4,k);  
    ab_snr1 = getsignalnoise(results.(whichmodel), 2, 'SNR'); % Left
    ab_snr2 = getsignalnoise(results.(whichmodel), 3, 'SNR'); % Right
    ab_snr_diff = to157chan(ab_snr2-ab_snr1,~data{1}.badChannels,'nans');
    
    [~,ch] = megPlotMap(ab_snr_diff,[-5,5],gcf,jmaColors('coolhotcortex'));
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-5:1:5);
    makeprettyaxes(ch,9,9);
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure9ab_bbRightMLeft_after'),[],0,'.',1);
end



