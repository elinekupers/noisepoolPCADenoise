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
    };

%conditionNamesAll = {'ON FULL','OFF FULL','ON LEFT','OFF LEFT','ON RIGHT','OFF RIGHT'};
conditionNamesAll = {'ON FULL','ON RIGHT','ON LEFT','OFF FULL','OFF RIGHT','OFF LEFT'};
conditionNames    = conditionNamesAll(conditionNumbers);

sessionDir      = sessionDirsAll{sessionnum};
rootDir = strrep(which('DFDrootpath.m'),'DFDrootpath.m','');
%megDataDir = fullfile(rootDir,'MEG','data');
megDataDir = fullfile(inputDataDir);