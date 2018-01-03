function dfdDenoiseWrapper(subjects, howToDenoise)
%
% dfdDenoiseWrapper(subjects, howToDenoise)
% ---------------------------------------------------------------- 
% INPUTS:
% -----------------
% subjects     : Vector of subject numbers one would like to denoise. See
%                   subject description in analysis/data for corresponding
%                   subject numbers
% howToDenoise : 1, 2, or 3, meaning:
%                   1) denoise with exactly 10 PC regressors
%                   2) denoise with each of 0 to 10 PC regressors
%                   3) various controls
%                   4) control that uses all channels in the noise pool
%                [default=1]
%
% DESCRIPTION: Wrapper function to denoise multiple MEG visual steady
% state data sets. The results of the denoising are written to file and can
% be then be used to reproduce several of the published figures from the
% paper the paper:
%
% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)

% Example 1: Denoising subject 1 from NYU dataset with exactly 10 PCs
%   dfdDenoiseWrapper(1, 1)
% Example 2: Denoising all subjects from NYU dataset with all PCs
%   dfdDenoiseWrapper(1:8, 2)
% Example 3: Denoising all subjects from CiNet dataset with 10 PCs
%   dfdDenoiseWrapper(9:12, 1)


% ------------------------------------------------------------------------
% --------------------------- Check options ------------------------------
% ------------------------------------------------------------------------

% Check for how to denoise the data
if notDefined('howToDenoise'), howToDenoise = 1; end

% Check for data, download data if needed
if isempty(dir(fullfile(dfdRootPath, 'analysis', 'data', 's*.mat')))
    error('No data were found. Use dfdDownloadSampleData')
end

% ------------------------------------------------------------------------
% --------------- Define variables depending on how to denoise -----------
% ------------------------------------------------------------------------

% Preprocessing parameters to remove bad channels and epochs (see dfdPreprocessData)
varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
opt.verbose         = true;
removeFirstEpoch    = true;  if removeFirstEpoch; opt.removeFirstEpoch = 1; end
removeMsEpochs      = false; if removeMsEpochs;   opt.removeMsEpochs   = 1; end %#ok<UNRCH>

%% Get frequencies to define stimulus locked and asynchronous broadband power
% Data are sampled at 1000 Hz and epoched in one-second bins. Hence
% frequencies are resolved in 1 Hz increments, 0:999 Hz
f           = 0:150;   % limit frequencies to [0 150] Hz
sl_freq     = 12;      % Stimulus-locked frequency
sl_freq_i   = sl_freq + 1;
tol         = 1.5;     % exclude frequencies within +/- tol of sl_freq
sl_drop     = f(mod(f, sl_freq) <= tol | mod(f, sl_freq) > sl_freq - tol);
   
% Exclude all frequencies that are close to a multiple of the
% line noise frequency (60 Hz)
ln_drop     = f(mod(f, 60) <= tol | mod(f, 60) > 60 - tol);

% Exclude all frequencies below 60 Hz when computing broadband power
lf_drop     = f(f<60);

% Define the frequenies and indices into the frequencies used to compute
% broadband power
[~, ab_i]   = setdiff(f, [sl_drop ln_drop lf_drop]);

keep_frequencies    = @(x) x(ab_i);
bb_frequencies      = f(ab_i);

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,sl_freq_i); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,keep_frequencies,1000);  % function handle to compuite broadband with a sample rate of 1 kHz

switch howToDenoise % Define denoise other parameters (see denoisedata.m)
    case {1, 4} % Denoise with exactly 10 PC regressors
        opt.pcchoose          = -10;   % Take 10 PCs
        opt.npcs2try          = [];
        optsl                 = opt;
        optbb                 = opt;
        optbb.preprocessfun   = @(x) bbFilter(x, bb_frequencies);  % preprocess data with a filter for broadband analysis
        nrControlModes        = 0;
        postFix               = '';
        saveDenoisedts        = true;
        
    case 2 % Denoise with each of 0 to 10 PC regressors
        opt.pcchoose          = 1.05;   % Get threshold for optimal nr of PCs
        opt.npcs2try          = '';     % empy string means up to nr of channels in noise pool
        optsl                 = opt;
        optbb                 = opt;
        optbb.preprocessfun   = @(x) bbFilter(x, bb_frequencies);  % preprocess data with a filter for broadband analysis
        nrControlModes        = 0;
        postFix               = '_full';
        saveDenoisedts        = false;
        
    case 3 % Denoise with various control modes
        opt.pcchoose          = -10;     % Take 10 PCs
        opt.npcs2try          = []; 
        optsl                 = opt;
        optbb                 = opt;
        optbb.preprocessfun   = @(x) bbFilter(x, bb_frequencies);  % preprocess data with a filter for broadband analysis
        nrControlModes        = 1:5;   % All control modes
        postFix               = 'control';  
        saveDenoisedts        = false;
        
end

% ------------------------------------------------------------------------
% ------------------ Load and denoise data per subject -------------------
% ------------------------------------------------------------------------

for whichSubject = subjects
    
    if whichSubject < 9,        dataChannels        = 1:157; % yokogawa MEG 157 channels
    elseif whichSubject < 21,   dataChannels        = 1:204; % Elekta Neuromag
    elseif whichSubject == 99;  dataChannels        = 1:157; % Synthetic subject
    elseif whichSubject > 36,   dataChannels        = 1:208; % yokogawa MEG 208 channels
    else                        dataChannels        = 1:157;  %#ok<SEPEX>
    end

    % ------------------ Load data and design ----------------------------
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    if removeMsEpochs
    tmp = load(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 'eye', 's%02d_epochswithms.mat'),whichSubject)); msepochidx = cat(1,tmp.allEpochsWithMS{:}); end %#ok<UNRCH>


    % ------------------ Make design matrix ------------------------------
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is left field
    design(conditions==7,3) = 1; % condition 7 is right field
    % condition 3 is blank
    
    % ------------- Combine channels if CiNet data -------------
    if length(dataChannels) == 157
        optbb.npoolmethod = {'r2','n',75};
        optsl.npoolmethod = {'r2','n',75};
    else
        optbb.npoolmethod = {'r2','n',100};
        optsl.npoolmethod = {'r2','n',100};
    end
    
    % ------------------ Preprocess data ---------------------------------
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold);
    
    % ---- Define first epochs in order to remove later ------------------
    if removeFirstEpoch, badEpochs(1:6:end) = 1; end
    
    % ---- Label epochs with microsaccades as bad epochs ------------------
    if removeMsEpochs 
        withMS = zeros(1,length(badEpochs)); %#ok<UNRCH>
        withMS(msepochidx)=1;
        badEpochs(find(withMS')) = 1; 
    end
     
    % -------------- Remove bad epochs and channels ----------------------
    sensorData      = sensorData(:,~badEpochs, ~badChannels);
    design          = design(~badEpochs,:);
    
    % ------------------ Permute sensorData for denoising ----------------
    sensorData = permute(sensorData, [3 1 2]);  
    
% ------------------------------------------------------------------------
% -------------------- Denoise and save the data -------------------------
% ------------------------------------------------------------------------
if howToDenoise == 4
    optbb.epochgroup = ones(1,size(sensorData,3));
    optsl.epochgroup = ones(1,size(sensorData,3));
    postFix          = '_denoise_all';   
end
    
    % ------------------ Denoise for broadband analysis ------------------
    for nrControl = nrControlModes   
        
        optbb.preprocessfun   = @(x) bbFilter(x, bb_frequencies);

        if (0 <= nrControl) && (nrControl <= 4)
            optbb.pccontrolmode   = nrControl;
            if nrControl == 0     % do nothing to postFix in filename
            else postFix          = sprintf('_control%d',nrControl); end %#ok<SEPEX>
            [results,evalout,~,denoisedts_bb] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);  %#ok<ASGLU>
            if saveDenoisedts; save(sprintf(fullfile(dfdRootPath, 'analysis', 'data', 's%02d_denoisedts.mat'),whichSubject),'denoisedts_bb'); end
            
        elseif nrControl == 5
            optbb.pccontrolmode   = 0;
            optbb.npoolmethod     = {'snr','n',size(sensorData,1)};
            postFix               = sprintf('_control%d',nrControl);
            [results,evalout] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
        end
  
        % ------------------ Define file name ---------------------------
        fname = sprintf(fullfile(dfdRootPath,'analysis','data', ['s%02d_denoisedData' postFix]), whichSubject);   
            
        % ----------------- Save denoised broadband data -----------------
        results.opt.preprocessfun   = func2str(results.opt.preprocessfun);
        optbb.preprocessfun         = func2str(optbb.preprocessfun);
        parsave([fname '_bb.mat'], 'results', results, 'evalout', evalout, 'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)
                
        % ----------- Denoise and save stimulus-locked analysis ----------
        [results,evalout] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);
        parsave([fname '_sl.mat'], 'results', results, 'evalout', evalout, 'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optsl)
    end
    if opt.verbose, sprintf('(%s) done denoising dataset subject s%02d with control number %d\n',mfilename,whichSubject,nrControl); end
    
    
end
