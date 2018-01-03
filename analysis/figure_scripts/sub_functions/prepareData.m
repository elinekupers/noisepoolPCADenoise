function [data,design,exampleIndex,exampleChannel] = prepareData(dataDir,whichSubject,whichFigure)

switch whichFigure
    case {4, 6}
        % Load denoised data
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_sensorData.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_denoisedts.mat'),whichSubject));
        
        exampleChannel = 42;
        
        % preprocessing parameters (see dfdPreprocessData)
        varThreshold        = [0.05 20];
        badChannelThreshold = 0.2;
        badEpochThreshold   = 0.2;
        if whichSubject < 9 || whichSubject == 99
            dataChannels    = 1:157;
        else
            dataChannels    = 1:204;
        end
        
        use3Channels        = false;
        
        % Preprocess raw sensordata
        [sensorData, badChannels0, badEpochs0] = dfdPreprocessData(sensorData(:,:,dataChannels), ...
            varThreshold, badChannelThreshold, badEpochThreshold, use3Channels);
        
        % ---- Define first epochs in order to remove later ------------------
        badEpochs0(1:6:end) = 1; 

        % Remove bad channels and bad epochs from data and conditions
        sensorData = sensorData(:,~badEpochs0, ~badChannels0);
                 
        % Permute sensorData for denoising
        sensorData = permute(sensorData, [3 1 2]);
        
        % time domain data before and after denoising
        denoisedts = denoisedts_bb;
        data = {sensorData,denoisedts{1}}; 
        
    case {5, 12, 14, 'SF5'}
        % Load denoised stimulus locked and broadband data
        bb = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject)); 
        sl = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_sl.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
        data = {bb,sl};
        
    case 7
        data = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_full_bb.mat'),whichSubject));        
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
        
    case {8, 9, 10, 'SF6'}
        data = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
        
    case 11
        for nrControl = 1:5
            data_controls{nrControl} = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_control%d_bb.mat'),whichSubject,nrControl)); %#ok<AGROW>
        end
        data_controls{6} = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_denoise_all_bb.mat'),whichSubject)); 
        
        data = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject)); 
        data = {data,data_controls};
        
   
    case 'SF1'
        data = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_NCPSvsNoisePool_bb.mat'), whichSubject));
        data1 = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
        data.badEpochs = data1.badEpochs;
        data.badChannels = data1.badChannels;

        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));
        
    case 'SF2'
        data = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_varyEpochLength_NrPCs_bb.mat'), whichSubject));
        data1 = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_bb.mat'),whichSubject));
        data.badEpochs = data1.badEpochs;
        data.badChannels = data1.badChannels;
        
        clear data1

        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));     
        
    case 'SF7'
        data = load(sprintf(fullfile(dataDir, 's%02d_denoisedData_sl.mat'),whichSubject));
        load(sprintf(fullfile(dataDir, 's%02d_conditions.mat'),whichSubject));     

end

%% Spectrum before and after denoising

% Because some channels were removed prior to the denoising (due to
% excessive noise), there are fewer channels in the denoised data than in
% the original data. We need to identify the correspondence.
if ~exist('exampleChannel','var')
    exampleIndex = [];
    
    
    % Compute design matrix
    design = zeros(size(conditions,1),3);
    design(conditions == 1,1) = 1; % Full
    design(conditions == 5,2) = 1; % Right
    design(conditions == 7,3) = 1; % Left
    
    if whichFigure == 4
        design = [];
%     elseif whichFigure == 14
%         design = [];
    elseif iscell(data)
        design = design(~data{1}.badEpochs,:);
    else
        design = design(~data.badEpochs,:);
    end
    
else
    exampleIndex = megGetOrigChannel(exampleChannel,badChannels);
    
    % Compute design matrix
    design = zeros(size(conditions,1),3);
    design(conditions == 1,1) = 1; % Full
    design(conditions == 5,2) = 1; % Right
    design(conditions == 7,3) = 1; % Left
    
    design0 = design(~badEpochs0,:);
    design1 = design(~badEpochs,:);
    design = {design0, design1};
end






