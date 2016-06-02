

%% s_SSMEG_TSPCA
%
% This script loads preprocesssed data from SSMEG experiment, conducts timeshifted PCA, and resaves.
% The saved files from this script have the names:
%   s0*_conditions.mat
%   s0*_sensorData.mat
% Where s0* refers to each of subjects 1-8

% Depends on Meg Utils
addpath(genpath('~/matlab/meg_utils'));

%% -------------------------------------
% -------------- Load data -------------
% --------------------------------------

project_pth                   = '/Volumes/server/Projects/MEG/Gamma_BR/data';
subjects                      = 2;
data_channels                 = 1:157;
trigger_channels              = 161:164;
photodiode_channel            = 192; 

fs                            = 1000;        % sample rate (Hz)
epoch_time                    = [0 1];       % start and end of epoch, relative to onset, in seconds
%% To run script, you need the Field trip toolbox

% Add fieldtrip path
field_trip_pth = fullfile(fileparts(project_pth), 'code', 'fieldtrip');
% meg_add_fieldtrip_paths(field_trip_pth, 'yokogawa_defaults');
addpath(genpath(field_trip_pth))

% Find subjects for this project
subject_pths = dir(fullfile(project_pth,'*Gamma*'));

% Make sure there are exactly 8 subjects
% assert(length(subject_pths)==8)

%% -------------------------------------
% ------- STRUCTURE THE RAW DATA -------
% --------------------------------------
%% Load in data
save_pth = fullfile(project_pth, 'manuscript_data');
% if ~exist(save_pth, 'dir'), mkdir(save_pth); end

for which_subject = subjects
    data_pth = fullfile(project_pth, subject_pths(which_subject).name, 'raw');
    
    subj_pth    = '*Gamma*';
    d = dir(fullfile(data_pth, subj_pth));
    
    rawFile   = fullfile(project_pth, subject_pths(which_subject).name, 'raw', d(1).name);
    CALMFile = fullfile(project_pth, subject_pths(which_subject).name, 'raw', d(2).name);  
    
%     data = meg_load_sqd_data(fullfile(project_pth, subject_pths(which_subject).name, 'raw'),'*SSMEG*');
    
%     sqdwrite(rawFile,tspcaFile,data);
%     tspcaFile1 = fullfile(project_pth, subject_pths(which_subject).name, 'raw', d(3).name); 
%     tspcaFile2 = fullfile(project_pth, subject_pths(which_subject).name, 'raw', d(3).name); 
    



%% -------------------------------------
% ------------ Perform TSPCA -----------
% --------------------------------------

%% Time-shift PCA for environmental denoising
% http://www.isr.umd.edu/Labs/CSSL/simonlab/Denoising.html
% http://lumiere.ens.fr/Audition/adc/meg/
% "Noise fields measured by reference magnetometers are optimally filtered 
% and subtracted from brain channels. The filters (one per reference/brain 
% sensor pair) are obtained by delaying the reference signals, 
% orthogonalizing them to obtain a basis, projecting the brain sensors onto 
% the noise-derived basis, and removing the projections to obtain clean 
% data." (DOI: 10.1016/j.jneumeth.2007.06.003)
    sizeOfBlocks = 20000;
    shifts = 1;%-100:100;

    sqdDenoise(sizeOfBlocks, shifts, 0, rawFile, [18 20  115  152], 'no', 168, 'no', CALMFile);
    
%
end
    






% %% -------------------------------------
% % --------------- Save Data ------------
% % --------------------------------------
% 
% fname = sprintf(fullfile(dfdRootPath,'exampleAnalysis','data', ['s%02d_sensorData_TSPCA_%d']),whichSubject, shifts);
% parsave([fname '.mat'], 'sensorData', sensorData)

