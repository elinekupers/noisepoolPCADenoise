function dfdPreprocessData()

%% Description etc



%% Variables
subjects                = 1:8;
inputDataDir            = fullfile(DFDrootpath, 'data', 'raw');
outputDataDir           = fullfile(DFDrootpath, 'data', 'preprocessed'); 
sensorDataPths          = dir(fullfile(inputDataDir, '*sensorData*'));
conditionsPths          = dir(fullfile(inputDataDir, '*conditions*'));
opt.remove_badepochs    = true;
data_channels           = 1:157;

data                    = [];

%% Loop over subjects
for whichSubject = subjects
    
    %% Load raw data
    try
        data{whichSubject,1} = load(fullfile(inputDataDir, sensorDataPths(whichSubject).name));
        data{whichSubject,2} = load(fullfile(inputDataDir, conditionsPths(whichSubject).name));
        
    catch ME
        idSegLast = regexp(ME.identifier, 'Undefined');
        if idSegLast
            fprintf('(dfdPreProcessData): ERROR. Data cannot be found. Please download data with the function dfdDownloadData');
        else
            disp(ME)
        end
    end
    
    %% Make design matrix
    condition_numbers = unique(data{whichSubject,2}.conditions);
%     assert(condition_numbers == [1;3;5;7])

    design = zeros(length(data{whichSubject,2}.conditions), 3);
    design(data{whichSubject,2}.conditions==1,1) = 1;
    design(data{whichSubject,2}.conditions==5,2) = 1;
    design(data{whichSubject,2}.conditions==7,3) = 1;
    
    if ~exist('outputDataDir','dir'), mkdir(outputDataDir); end
    save(fullfile(outputDataDir, sprintf('s%02d_design', whichSubject)), 'design');
    
    %% Remove bad channels
    sensorData = data{whichSubject,1}.sensorData(:,:,1:157);
    
    opt.remove_badepochs = true;
    
        okEpochs = true(size(sensorData,2),1);
        [sensorData, okEpochs] = dfdIdenitfyBadEpochs(sensorData, data_channels, 0.5);
%         okEpochs = okEpochs & ;
        
        sensorData = sensorData(:,:,okEpochs);
        design = design(okEpochs,:);
    
    
    
end











%% Remove bad epochs



%% Save prepocessed data



