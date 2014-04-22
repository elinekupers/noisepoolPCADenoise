function [sessionDir, conditionNames, megDataDir] = megGetDataPaths(sessionum, conditionNumbers)

sessionDirsAll = {...
    '02_SSMEG_02_28_2014';...
    '03_SSMEG_03_31_2014';...
    '04_SSMEG_04_01_2014';...
    '05_SSMEG_04_04_2014'
    };

conditionNamesAll = {'ON FULL','OFF FULL','ON LEFT','OFF LEFT','ON RIGHT','OFF RIGHT'};

sessionDir      = sessionDirsAll{sessionum};
conditionNames  = conditionNamesAll(conditionNumbers);

rootDir = strrep(which('setup.m'),'denoisesuite/setup.m','');
megDataDir = fullfile(rootDir,'MEG','data');