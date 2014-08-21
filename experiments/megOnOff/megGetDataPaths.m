function [sessionDir, megDataDir,conditionNames] = megGetDataPaths(sessionum,conditionNumbers)
if notDefined('conditionNumbers'), conditionNumbers = 1:6; end

sessionDirsAll = {...
    '02_SSMEG_02_28_2014_wl_subj004'; ... 1
    '03_SSMEG_03_31_2014_wl_subj005'; ... 2
    '04_SSMEG_04_01_2014_wl_subj006'; ... 3
    '05_SSMEG_04_04_2014_wl_subj002'; ... 4
    '06_SSMEG_04_28_2014_wl_subj008'; ... 5
    '07_SSMEG_05_01_2014_wl_subj009'; ... 6
    '08_SSMEG_06_20_2014_wl_subj011'; ... 7
    '09_SSMEG_06_27_2014_wl_subj010'; ... 8
    ...
    '05_SSMEG_Attention_07_16_2014_wl_subj002'; ... 9
    '06_SSMEG_Attention_07_17_2014_wl_subj010'; ... 10
    '07_SSMEG_Attention_07_21_2014_wl_subj012'; ... 11
    '08_SSMEG_Attention_07_23_2014_wl_subj005'; ... 12
    '09_SSMEG_Attention_07_31_2014_wl_subj011'; ... 13
    '10_SSMEG_Attention_08_04_2014_wl_subj014'; ... 14
    ...
    };

%conditionNamesAll = {'ON FULL','OFF FULL','ON LEFT','OFF LEFT','ON RIGHT','OFF RIGHT'};
conditionNamesAll = {'ON FULL','ON RIGHT','ON LEFT','OFF FULL','OFF RIGHT','OFF LEFT'};
conditionNames    = conditionNamesAll(conditionNumbers);

sessionDir      = sessionDirsAll{sessionum};
rootDir = strrep(which('setup.m'),'denoisesuite/setup.m','');

% megDataDir = fullfile(rootDir,'MEG','data');
megDataDir = '/Volumes/server/Projects/MEG/';
if strfind(sessionDir, 'Attention')
    megDataDir = fullfile(megDataDir, 'Attention_MEG', 'Data');
else
    megDataDir = fullfile(megDataDir, 'SSMEG');
end