function dfdPreprocessData(subjects)

%% Description etc



%% Variables
if notDefined('subjects'), subjects = 1:8; end
inputDataDir            = fullfile(dfdRootPath, 'data', 'raw');
outputDataDir           = fullfile(dfdRootPath, 'data', 'preprocessed'); 
sensorDataPths          = dir(fullfile(inputDataDir, '*sensorData*'));
conditionsPths          = dir(fullfile(inputDataDir, '*conditions*'));
opt.removeBadepochs    = true;
dataChannels           = 1:157;



%% Loop over subjects

for whichSubject = subjects
    
    %% Load raw data
    sensorDataPth = fullfile(inputDataDir, sensorDataPths(whichSubject).name);
    conditionsPth = fullfile(inputDataDir, conditionsPths(whichSubject).name);
    if exist(sensorDataPth, 'file') && exist(conditions_pth, 'file')
        load(sensorDataPth); % assumption: sensorDataPth is a file with a variable called 'sensorData'
        load(conditionsPth); % assumption: conditionsPth is a file with a variable called 'conditions'
    else        
        fprintf('(dfdPreProcessData): ERROR. Data cannot be found. Please download data with the function dfdDownloadData');
    end
    
    %% Make design matrix
    condition_numbers = unique(conditions);
%     assert(condition_numbers == [1;3;5;7])

    glmDesign = zeros(length(conditions), 3);
    glmDesign(conditions==condition_numbers(1),1) = 1;
    glmDesign(conditions==condition_numbers(3),2) = 1;
    glmDesign(conditions==condition_numbers(4),3) = 1;
    
    if ~exist('outputDataDir','dir'), mkdir(outputDataDir); end
    save(fullfile(outputDataDir, sprintf('s%02d_design', whichSubject)), 'design');
    
    %% Remove bad channels
    sensorData = sensorData(:,:,dataChannels);
    
    
    okEpochs = true(size(sensorData,2),1);
    [sensorData, okEpochs] = dfdIdenitfyBadEpochs(sensorData, dataChannels, 0.5);
    %         okEpochs = okEpochs & ;
    
    sensorData = sensorData(:,:,okEpochs);
    glmDesign = glmDesign(okEpochs,:);
    
    
    
end











%% Remove bad epochs



%% Save prepocessed data



