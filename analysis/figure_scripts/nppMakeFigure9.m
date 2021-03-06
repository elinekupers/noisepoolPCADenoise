function fH = nppMakeFigure9(whichSubject)
%% Function to reproduce Figure 9 (Spatialmap) from example subject and across subjects
% after denoising
%
% nppMakeFigure9()
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (2018) A non-invasive, quantitative study of
% broadband spectral responses in human visual cortex (PLOS ONE. VOLUME.
% ISSUE. DOI.)
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the broadband signals before and
% after using the denoising algorithm. The four separate conditions (Full,
% left, right and left-right hemifield stimulation are shown separately). 
%
% This function assumes that data is downloaded with the nppDownloaddata
% function. 

%% Choices to make: 

if ~exist('whichSubject', 'var')
    whichSubject    = 1;        % Subject 1 is the example subject.
end
figureDir       = fullfile(nppRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(nppRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;        % Set threshold for colormap. If no threshold set value to 0
cfg             = [];
data_hdr        = [];


%% Load denoised data of example subject
bb = prepareData(dataDir,whichSubject,9);

%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
condNames = {'Stim Full','Stim Left','Stim Right' 'Stim Left - Right'};
contrasts = [eye(3); [0 1 -1]/sqrt(2)];
yscaleAB = [repmat([-8,-4,0,4,8],3,1);[-5,-2.5,0,2.5,5]];

fH = figure('position',[1,600,1400,800], 'Name', 'Figure 9, Subject 1', 'NumberTitle', 'off');
for icond = 1:4
    
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bb.results.origmodel(1),  contrasts(icond,:), 'SNR',bb.badChannels);
    ab_snr2 = getsignalnoise(bb.results.finalmodel(1), contrasts(icond,:), 'SNR',bb.badChannels);
    

    % Replace bad channels with NaNs
    ab_snr1 = to157chan(ab_snr1,~bb.badChannels,'nans');
    ab_snr2 = to157chan(ab_snr2,~bb.badChannels,'nans');
                
    if length(bb.badChannels)>157 && length(bb.badChannels) <= 204
        ab_snr1 = npp204to102(ab_snr1);
        ab_snr2 = npp204to102(ab_snr2);
    end        
    
    % Threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    ab_snr2(abs(ab_snr2) < threshold) = 0;
    
    
    % Set colormap limits
%     max_val = max(abs([ab_snr1a_LmnR, ab_snr2a_LmnR]));
%     clims_ab = [-1,1].*[max_val,max_val];
    if icond <= 3, clims_ab = 8.4445 * [-1 1]; else, clims_ab = 3.5363 * [-1,1]; end
    
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(ab_snr1,clims_ab,gcf,'bipolar',sprintf('%s Original', condNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',yscaleAB(icond,:));
    title(sprintf('Broadband Pre %s', condNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr2,clims_ab,gcf,'bipolar',sprintf('%s : Denoised PC %d',condNames{icond}, bb.results.pcnum(1)),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',yscaleAB(icond,:));
    title(sprintf('Broadband Post %s', condNames{icond}))
    drawnow();
    
end

if saveFigures
      % Only use figure write when producing MS high quality figure, since
      % it takes time
     hgexport(gcf,fullfile(figureDir, sprintf('figure9_examplesubject%d_bipolar_thresh%d',whichSubject, threshold)));
%     figurewrite(fullfile(figureDir, sprintf('figure9_examplesubject%d_bipolar_thresh%d_diff',whichSubject, threshold)),[],0,'.',1);
end

%% Now do the same but then across subjects
if whichSubject ~= 99
    nppMakeFigure9AcrossSubjects()
end
