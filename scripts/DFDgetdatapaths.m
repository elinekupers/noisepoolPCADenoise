function [sessionDir, megDataDir,conditionNames] = DFDgetdatapaths(sessionnum,conditionNumbers, inputDataDir)
if notDefined('conditionNumbers'), conditionNumbers = 1:6; end

sessionDirsAll = {...
    'SSMEG_Dataset_04'; %1
    'SSMEG_Dataset_05'; %2
    'SSMEG_Dataset_06'; %3
    'SSMEG_Dataset_07'; %4
    'SSMEG_Dataset_08'; %5
    'SSMEG_Dataset_09'; %6
    'SSMEG_Dataset_10'; %7
    'SSMEG_Dataset_11'; %8

% %     '02_SSMEG_02_28_2014';... %1
% %     '03_SSMEG_03_31_2014';... %2
%     '04_SSMEG_04_01_2014';... %3
%     '05_SSMEG_04_04_2014';... %4
%     '06_SSMEG_04_28_2014';... %5
%     '07_SSMEG_05_01_2014';... %6
%     ...
% %     '01_SSMEG_Attention_wl_subj002'; ... %7
% %     '02_SSMEG_Attention_wl_subj010'; ... %8
%     ...
%     '08_SSMEG_06_20_2014_subj011';... %9,
%     '09_SSMEG_06_27_2014_subj010';... %10
%     ...
%     '10_SSMEG_08_12_2014_wl_subj004';...%this replaces 1
%     '11_SSMEG_08_13_2014_wl_subj005';...%this replaces 2
    
    };

%conditionNamesAll = {'ON FULL','OFF FULL','ON LEFT','OFF LEFT','ON RIGHT','OFF RIGHT'};
conditionNamesAll = {'ON FULL','ON RIGHT','ON LEFT','OFF FULL','OFF RIGHT','OFF LEFT'};
conditionNames    = conditionNamesAll(conditionNumbers);

sessionDir      = sessionDirsAll{sessionnum};
rootDir = strrep(which('DFDrootpath.m'),'DFDrootpath.m','');
%megDataDir = fullfile(rootDir,'MEG','data');
megDataDir = fullfile(inputDataDir);