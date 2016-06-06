function dfdMakeFigure8()
%% Function to reproduce Figure 8 (Spatialmap) for all subjects pre- vs post-denoising
%
% dfdMakeFigure8()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in
% each channel for the broadband signals before and
% after using the denoising algorithm. Only full field condition is shown.
%
% This function assumes that data is downloaded with the dfdDownloadData
% function.

%% Choices to make:

whichSubjects        = 1:8;
dataDir              = fullfile(dfdRootPath, 'exampleAnalysis','data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'exampleAnalysis','figures');% Where to save images?
saveFigures          = true;   % Save figures in the figure folder?
dataAll              = [];

%% Load data of all subjects
for whichSubject = whichSubjects
    fprintf(' Load subject %d \n', whichSubject);
    [data,design,exampleIndex] = prepareData(dataDir,whichSubject,8);
    dataAll{whichSubject} = {data,design,exampleIndex}; %#ok<AGROW>
end

%% Plot spatial map figure with full condition, post minus pre denoising
figure('position',[1,600,1400,800]);
icond = 1;
for k = 1:length(whichSubjects)
    
    data = dataAll{k};
    results = data{1}.results;
    
    subplot(2,4,k);
    ab_snr1 = getsignalnoise(results.origmodel,  icond, 'SNR');
    ab_snr2 = getsignalnoise(results.finalmodel, icond, 'SNR');
    ab_snr_diff = to157chan(ab_snr2-ab_snr1,~data{1}.badChannels,'nans');
    
    [~,ch] = megPlotMap(ab_snr_diff,[-10,10],gcf,jmaColors('coolhotcortex'));
    makeprettyaxes(gca,9,9);
    set(ch,'ytick',-10:5:10);
    makeprettyaxes(ch,9,9);
end

if saveFigures
    figurewrite(fullfile(figureDir,'figure8_bbdiffall'),[],0,'.',1);
end
return

