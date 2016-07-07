function s_ERF_analysis()

% Script to plot the eye related function of the data

subjects = 8;

% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];%[0.05 30];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
use3Channels        = true;
removeFirstEpoch    = true;
triggerChannels    = 161:164;
environmentalChannels = 158:160;
verbose             = true;
% 
% 
% %% Get frequencies to define stimulus locked and asynchronous broadband power
% % Data are sampled at 1000 Hz and epoched in one-second bins. Hence
% % frequencies are resolved in 1 Hz increments, 0:999 Hz
% f           = 0:150;   % limit frequencies to [0 150] Hz
% sl_freq     = 12;      % Stimulus-locked frequency
% sl_freq_i   = sl_freq + 1;
% tol         = 1.5;     % exclude frequencies within +/- tol of sl_freq
% sl_drop     = f(mod(f, sl_freq) <= tol | mod(f, sl_freq) > sl_freq - tol);
%    
% % Exclude all frequencies that are close to a multiple of the
% % line noise frequency (60 Hz)
% ln_drop     = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);
% 
% % Exclude all frequencies below 60 Hz when computing broadband power
% lf_drop     = f(f<60);
% 
% % Define the frequenies and indices into the frequencies used to compute
% % broadband power
% [~, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);
% 
% keep_frequencies    = @(x) x(ab_i);
% bb_frequencies      = f(ab_i);
% 
% % Define functions to define noise pool and signal of interest
% evokedfun           = @(x)getstimlocked(x,sl_freq_i); % function handle to determine noise pool
% evalfun             = @(x)getbroadband(x,keep_frequencies,1000);  % function handle to compuite broadband with a sample rate of 1 kHz
% 
% % Define options for denoising that are equal for each type of denoising
% opt.resampling      = {'boot','boot'};
% opt.pcselmethod     = 'snr';
% opt.verbose         = true;

project_pth                   = '/Volumes/server/Projects/MEG/SSMEG';

%% To run script, you need the Field trip toolbox
% 
% % Add fieldtrip path
% field_trip_pth = fullfile(fileparts(project_pth), 'code', 'fieldtrip');
% % meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');
% addpath(genpath(field_trip_pth))

% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*SSMEG_*'));


for whichSubject = subjects
    
    if whichSubject < 9,        dataChannels        = 1:157; % yokogawa MEG
    elseif whichSubject < 21,   dataChannels        = 1:204; % Elekta Neuromag
    else                        dataChannels        = 1:157; % Synthetic subject
    end

    % ------------------ Load ms timings ----------------------------
%     tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 'eye', 's%02d_ms.mat'),whichSubject)); ms = tmp.ms;
    
%     % ------------------ Load raw ts ----------------------------
%     data_pth = fullfile(project_pth, subject_pths(whichSubject).name, 'raw');
%     [tsRaw, meg_files] = meg_load_sqd_data(data_pth, '*SSMEG_*');
    

    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_sensorDataERF.mat'),whichSubject)); tsRaw = tmp.sensorData;
    
    
    % ------------------ Load denoised ts ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_denoisedts_ERF.mat'),whichSubject)); tsDenoised = tmp.denoisedts_bb{1};
    tsDenoised = permute(tsDenoised,[2,3,1]);
    
    tsRaw = dfdEnvironmentalDenoising(tsRaw, environmentalChannels, dataChannels);

    
%     figure;
%     for chan = 1:157
%         clf;
%         subplot(211)
%         plot(mean(tsRaw(:,:,chan),2),'b'); hold on;
%         plot(mean(tsDenoised(:,:,chan),2),'k');
%         xlabel('Time from ms onset [ms]')
%         ylabel('MS locked ERF (fTesla)')
%         legend('Raw','Denoised')
%        
%         subplot(212)
%         electrodes = zeros(1,157);
%         electrodes(chan)=1;
%         ft_plotOnMesh(electrodes, chan, [], [],'interpolation','nearest')
%         
%         pause(1);
%     end
%         
    fftRaw = abs(fft(mean(tsRaw(1:200,:,chan),2)));
        
    figure; hold on;
    
    subplot(811); 
    for chan = [9,11,12,79,75,78]   
        hold all;
        plot(mean(tsRaw(1:200,:,chan),2),'b');
        plot(mean(tsDenoised(1:200,:,chan),2),'k');
        xlabel('Time from ms onset [ms]')
        ylabel('MS locked ERF (fTesla)')
        legend('Raw','Denoised')
    end
    
    subplot(812); 
    electrodes = zeros(1,157);
    electrodes([9,11,12,79,75,78])=1;
    ft_plotOnMesh(electrodes, chan, [], [],'interpolation','nearest')
    
    subplot(813);
    for chan = [49,40,38,56,53,52,37]
         hold all;
        plot(mean(tsRaw(1:200,:,chan),2),'b');
        plot(mean(tsDenoised(1:200,:,chan),2),'k');
        xlabel('Time from ms onset [ms]')
        ylabel('MS locked ERF (fTesla)')
        legend('Raw','Denoised')
    end
    subplot(814);
    electrodes = zeros(1,157);
    electrodes([49,40,38,56,53,52,37])=1;
    ft_plotOnMesh(electrodes, chan, [], [],'interpolation','nearest')
    

    subplot(815); 
    for chan = [108,109,155,122,123]
        hold all;
        plot(mean(tsRaw(1:200,:,chan),2),'b');
        plot(mean(tsDenoised(1:200,:,chan),2),'k');
        xlabel('Time from ms onset [ms]')
        ylabel('MS locked ERF (fTesla)')
        legend('Raw','Denoised')
    end
    subplot(816);
    electrodes = zeros(1,157);
    electrodes([108,109,155,122,123])=1;
    ft_plotOnMesh(electrodes, chan, [], [],'interpolation','nearest')
    
    
    subplot(817); 
    for chan = [117,119,106,102, 115];
        hold all;
        plot(mean(tsRaw(1:200,:,chan),2),'b');
        plot(mean(tsDenoised(1:200,:,chan),2),'k');
        xlabel('Time from ms onset [ms]')
        ylabel('MS locked ERF (fTesla)')
        legend('Raw','Denoised')
    end
    subplot(818);
    electrodes = zeros(1,157);
    electrodes([117,119,102,115,134])=1;
    ft_plotOnMesh(electrodes, chan, [], [],'interpolation','nearest')
    
    
%     
%   
%     onsets = [];
%     for ii = 1:numel(ms)
%         tmp = [find(ms{ii}) ones(length(find(ms{ii})),1)*ii];
%         onsets = cat(1,[onsets;tmp]);
%     end
%     
%     [eyetsRaw, ~]   = meg_make_epochs(tsRaw, onsets(:,1), [0 .999], 1000, 'eye');
%     
%     clear tsRaw
% 
%     % ------------------ Preprocess data ---------------------------------
%     [sensorData, badChannels, badEpochs] = dfdPreprocessData(eyetsRaw(:,:,dataChannels), ...
%         varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
%     
%     
%     
%     
    
    
    % Plot ERFs    
    figure;
    for chan = 1:157
        clf;
    plot(squeeze(mean(eyetsRaw(:,onsets(:,2)==1,chan),2)))
    pause(1);
    end
    
    

    
    
%     % ---- Define first epochs in order to remove later ------------------
%     if removeFirstEpoch, badEpochs(1:6:end) = 1; end
%     
%     % -------------- Remove bad epochs and channels ----------------------
%     sensorData      = sensorData(:,~badEpochs, ~badChannels);
%     design          = design(~badEpochs,:);
%     
%     % ------------------ Permute sensorData for denoising ----------------
%     sensorData = permute(sensorData, [3 1 2]); 
    
end

return

function sensorDataDenoised = dfdEnvironmentalDenoising(sensorData,environmentalChannels, dataChannels)
% Make empty arrays for regressed 'clean' data
sensorDataDenoised = sensorData;

warning off stats:regress:RankDefDesignMat

% Start regression, keep residuals
for channel = dataChannels; 

        fprintf('[%s]: Channel %d\n', mfilename, channel);

    for epoch = 1:size(sensorData,2);
        [~,~,R] = regress(sensorData(:,epoch,channel),[squeeze(sensorData(:,epoch,environmentalChannels)) ones(size(sensorData,1),1)]);
        sensorDataDenoised(:,epoch,channel) = R;
        clear R 
    end
end
warning on stats:regress:RankDefDesignMat

return
