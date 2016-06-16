%% Check how many channels in the noise pool is best for MEG

function dfdFindOptimalNoisepool(subjects)

% Description: Wrapper function to denoise multiple MEG visual steady
% state data sets. The results of the denoising are written to file and can
% be then be used to reproduce several of the published figures from the
% paper the paper:
%
% AUTHORS. YEAR. TITLE. JOURNAL.

% Check for data, download data if needed
if isempty(fullfile(dfdRootPath, 'data'));
    error('No data were found. Use dfdDownloadSampleData')
end

% preprocessing parameters (see dfdPreprocessData)
varThreshold        = [0.05 20];  
badChannelThreshold = 0.2;       
badEpochThreshold   = 0.2;
dataChannels        = 1:157;
use3Channels        = false;

% denoise parameters (see denoisedata.m)
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

% Define functions to define noise pool and signal of interest
evokedfun           = @(x)getstimlocked(x,sl_freq_i); % function handle to determine noise pool
evalfun             = @(x)getbroadband(x,keep_frequencies,1000);


optbb.preprocessfun   = @bbFilter;  % preprocess data with a filter for broadband analysis
npcs                  = [5,10:10:70];
optbb.resampling        = {'boot','boot'};
optbb.pcselmethod       = 'snr';
optbb.nboot             = 100;


% Load and denoise data, one subject at a time
for whichSubject = subjects
    allResults = [];
    % ******* Load data and design *****************
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis','data', 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis','data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    % ******** Make design matrix *****************
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is right (??)
    design(conditions==7,3) = 1; % condition 7 is left (??)
    % condition 3 is blank
    
    % ******* Preprocess data **********************
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
    
    % Remove bad channels and bad epochs from data and conditions
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design = design(~badEpochs,:);

    % Permute sensorData for denoising
    sensorData = permute(sensorData, [3 1 2]);
    
    % ********* Denoise the data ********************   
    for jj = 1:length(npcs) % iterate through number of PCs projected out
        optbb.pcchoose = -npcs(jj);
        fprintf('\tnpcs = %d\n', npcs(jj));
        
        if jj == 1
            [results] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
            noisepooldef = results.noisepool;
        else
            [results] = denoisedata(design,sensorData,noisepooldef,evalfun,optbb);
        end
        allResults{whichSubject,jj} = results;
    end
        
%      
%     
%     %   Denoise for broadband analysis    
% %     [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
%     
%     if use3Channels
%         fname = sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData_w3chan'),whichSubject);
%     else
%         fname = sprintf(fullfile(dfdRootPath,'data','s0%d_denoisedData'),whichSubject);
%     end
    fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_noisepoolresults'),whichSubject);
    parsave([fname '_bb.mat'], 'results', results)        
    
%     %   Denoise for stimulus-locked analysis
%     [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);
%     parsave([fname '_sl.mat'], 'results', results, 'evalout', evalout, ...
%         'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
%         'badChannels', badChannels, 'badEpochs', badEpochs,  'opt', optsl)        
%     
end