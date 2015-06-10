function dfdDenoiseWrapper(subjects)
%
% dfdDenoiseWrapper(subjects)
%
% INPUTS:
% subjects  : Number of subjects one would like to denoise
%
% DESCRIPTION: Wrapper function to denoise multiple MEG visual steady
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

% Get 'freq' struct to define stimulus locked and broadband frequencies
%  This struct is needed as input args for getstimlocked and getbroadband
freq = megGetSLandABfrequencies(0:150, 1, 12);

% denoise parameters (see denoisedata.m)
optsl.pcchoose        = -10;   % denoise with exactly 10 PCs for stimulus locked
optbb.pcchoose        = -10;   % denoise with exactly 10 PCs for broadband
optsl.npcs2try        = 0;   % loop through 10
optbb.npcs2try        = 0;
optsl.resampling      = {'boot','boot'};
optbb.resampling      = {'boot','boot'};
optsl.pcselmethod     = 'snr';
optbb.pcselmethod     = 'snr';
optbb.preprocessfun   = @hpf;  % preprocess data with a high pass filter for broadband analysis
evokedfun             = @(x)getstimlocked(x,freq); % function handle to determine noise pool
evalfun               = @(x)getbroadband(x,freq);  % function handle to compuite broadband
nrControlModes        = 0; % Note, mode 1 has a bug.
allInNoisepool        = false;
varyEpochLength       = true; if varyEpochLength, 
                                epochDurs = [1,3,6,12,24,36,72,inf];
                                npcs      = [5,10:10:70]; end
                            


% Load and denoise data, one subject at a time
for whichSubject = subjects
    % ******* Load data and design *****************
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 's0%d_sensorData.mat'),whichSubject)); sensorData = tmp.sensorData;
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    % ******** Make design matrix *****************
    design = zeros(length(conditions), 3);
    design(conditions==1,1) = 1; % condition 1 is full field
    design(conditions==5,2) = 1; % condition 5 is left field
    design(conditions==7,3) = 1; % condition 7 is right field
    % condition 3 is blank
    
    % ******* Preprocess data **********************
    [sensorData, badChannels, badEpochs] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
        varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
    
    % Remove bad channels and bad epochs from data and conditions
    sensorData = sensorData(:,~badEpochs, ~badChannels);
    design = design(~badEpochs,:);
    
    % Permute sensorData for denoising
    sensorData = permute(sensorData, [3 1 2]);
    
    % Get number of channels if denoising with all channels in noise pool
    if allInNoisepool
        optbb.npoolmethod = {'r2','n',size(sensorData,1)};
    end
    
    % ********* Denoise the data ********************
    
    
    for nrControl = nrControlModes;
        optbb.pccontrolmode = nrControl; % Zero is no control method
        
        if varyEpochLength
            % Vary epoch lengths
            if isinf(epochDurs(ii))
                epochGroup = ones(size(sensorData,1),1);
            else
                epochGroup = megEpochGroup(okEpochs,epochDurs(ii),0);
            end
            optbb.epochGroup = epochGroup;
            fprintf('epochDur = %d; %d epochs\n', epochDurs(ii), max(epochGroup));
        
                clear results;
        
%             % Vary npcs
%             for jj = 1:length(npcs) % iterate through number of PCs projected out 
%                 opt.pcstop = -npcs(jj);
%                 fprintf('\tnpcs = %d\n', npcs(jj));
% 
%                 if jj == 1
%                     [results] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
%                     noisepooldef = results.noisepool;
%                 else
%                     [results] = denoisedata(design,sensorData,noisepooldef,evalfun,opt);
%                 end
%                 allResults{ii,jj} = results;
%             end
%         end
        
        
        % Denoise for broadband analysis
        [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,optbb);
        
        
        % Get file name
        if optbb.pcchoose == 0
            if use3Channels
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_w3chan'),whichSubject);
            else
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData'),whichSubject);
            end
        elseif optbb.pccontrolmode > 0
            if use3Channels
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_w3chan_control%d'),whichSubject,optbb.pccontrolmode);
            else
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_control%d'),whichSubject,optbb.pccontrolmode);
            end
        elseif allInNoisepool
            if use3Channels
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_w3chan_allinnp'),whichSubject);
            else
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_allinnp'),whichSubject);
            end
        else
            if use3Channels
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_w3chan_full'),whichSubject);
            else
                fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data','s0%d_denoisedData_full'),whichSubject);
            end
            
            
            
        end
        
        % ********* Save denoised broadband data ********************
        parsave([fname '_bb.mat'], 'results', results, 'evalout', evalout, ...
            'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
            'badChannels', badChannels, 'badEpochs', badEpochs, 'opt', optbb)
        
        % ********* Denoise and save stimulus-locked analysis *****************
        [results,evalout,denoisedspec,denoisedts] = denoisedata(design,sensorData,evokedfun,evokedfun,optsl);
        parsave([fname '_sl.mat'], 'results', results, 'evalout', evalout, ...
            'denoisedspec', denoisedspec, 'denoisedts', denoisedts,...
            'badChannels', badChannels, 'badEpochs', badEpochs,  'opt', optsl)
    end
end
