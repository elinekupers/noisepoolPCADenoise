function dfdMakeFigure5()
%% Function to reproduce Figure 5 (Spatialmap) from example subject 
%
% dfdMakeFigure5()
%
% AUTHORS. TITLE. JOURNAL. YEAR.
%
% This figure will show an interpolated spatial map of the SNR values in 
% each channel for the stimulus locked signal, broadband signals before and
% after using the denoising algorithm. The three separate conditions (Full,
% left, right hemifield stimulation are shown separately). 
%
% This function assumes that data is downloaded with the DFDdownloaddata
% function. 

%% Choices to make:                                              
whichSubject    = 1;        % Subject 1 is the example subject.
figureDir       = fullfile(dfdRootPath, 'analysis', 'figures'); % Where to save images?
dataDir         = fullfile(dfdRootPath, 'analysis', 'data');    % Where to save data?
saveFigures     = true;     % Save figures in the figure folder?
threshold       = 0;        % Set threshold for colormap. If no threshold set value to 0
cfg             = [];
data_hdr        = [];


%% Load denoised data of example subject
[data] = prepareData(dataDir,whichSubject,5);
bb = data{1};
sl = data{2};

%% Plot stimulus-locked signal, broadband before and after denoising on sensormap
figure('position',[1,600,1400,800]); set(gcf, 'Name', 'Figure 5, Example subject', 'NumberTitle', 'off');
condNames = {'Stim Full','Stim Left','Stim Right'};
contrastNames = [condNames 'Left minus Right'];
contrasts = [eye(3); 0 1 -1];
contrasts = bsxfun(@rdivide, contrasts, sqrt(sum(contrasts.^2,2)));
yscaleAB = [repmat([-8,-4,0,4,8],3,1);[-5,-2.5,0,2.5,5]];

for icond = 1:size(contrasts,1)
    % get stimulus-locked snr
    sl_snr1 = getsignalnoise(sl.results.origmodel(1),contrasts(icond,:), 'SNR',sl.badChannels);
    % get broadband snr for before and after denoising
    ab_snr1 = getsignalnoise(bb.results.origmodel(1),  contrasts(icond,:), 'SNR',bb.badChannels);
    


    if whichSubject < 9; % NeuroMag360 data is already converted when combining the channels
        % convert back into 157-channel space  
        sl_snr1 = to157chan(sl_snr1,~sl.badChannels,'nans');
        ab_snr1 = to157chan(ab_snr1,~bb.badChannels,'nans');
           
    end
    
    % Threshold
    ab_snr1(abs(ab_snr1) < threshold) = 0;
    sl_snr1(abs(sl_snr1) < threshold) = 0;
    
    
   %CHECK IF THIS IS THE SAME 
%     ab_beta2 = results.whichmodel.beta(1,:,:); %right
%     ab_diff = ab_beta2 - ab_beta1;
%     
%     diff_med = nanmedian(squeeze(ab_diff),2);
%     diff_se = nanstd(squeeze(ab_diff),[],2);
%         
%     ab_snr_diff = to157chan((diff_med./diff_se)',~data{1}.badChannels,'nans');
%     
%     subplot(2,4,k);  

    % Set colormap limits
%     max_val = max(abs([ab_snr1a_LmnR, ab_snr2a_LmnR]));
%     clims_ab = [-1,1].*[max_val,max_val];
    clims_sl = [-25.6723,25.6723];
    clims_ab = [-8.4445, 8.4445];
    if icond > 3 % then we are plotting l-r rather than one condition
        clims_ab = [-5.5363, 5.5363]; 
    end

    % plot spatial maps
    subplot(4,2,(icond-1)*2+1)
    [~,ch] = megPlotMap(sl_snr1,clims_sl,gcf,'bipolar',sprintf('%s : Stimulus Locked Original', contrastNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',[-20,-10,0,10,20]);
    title(sprintf('SL no DN %s', contrastNames{icond}))
    
    subplot(4,2,(icond-1)*2+2)
    [~,ch] = megPlotMap(ab_snr1,clims_ab,gcf,'bipolar',sprintf('%s Original', contrastNames{icond}),data_hdr,cfg);
    makeprettyaxes(ch,9,9);
    set(ch,'YTick',yscaleAB(icond,:));
    title(sprintf('Broadband Pre %s', contrastNames{icond}))
    

end

if saveFigures
    printnice(gcf, 0, figureDir, sprintf('figure5_examplesubject%d_bipolar_thresh%d_interpolated_noMS',whichSubject, threshold));
end

%% Now call dfdMakeFigure5AcrossSubjects
dfdMakeFigure5AcrossSubjects();


