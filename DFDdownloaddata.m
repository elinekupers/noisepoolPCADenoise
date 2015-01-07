function DFDdownloaddata(savePth)
% Download sample MEG data sets to be denoised by the 'Denoise Field Data'
% algorithm for the paper:
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% DFDdownloaddata
%
% Example: 
%   DFDdownloaddata('/Users/jon/tmp/')

% Argument check
if nargin == 0, savePth = pwd; end

% Site to retrieve the data
dirProject = 'http://psych.nyu.edu/winawerlab/denoiseFieldData';

% These are the subject data sets
dirSubjects = {...
    %'01_SSMEG_Attention_wl_subj002'...
    %'02_SSMEG_02_28_2014'...
    %'02_SSMEG_Attention_wl_subj010'...
    %'03_SSMEG_03_31_2014'...
    '04_SSMEG_04_01_2014'...
    '05_SSMEG_04_04_2014'...
    '06_SSMEG_04_28_2014'...
    '07_SSMEG_05_01_2014'...
    '08_SSMEG_06_20_2014_subj011'...
    '09_SSMEG_06_27_2014_subj010'...
    '10_SSMEG_08_12_2014_wl_subj004'...
    '11_SSMEG_08_13_2014_wl_subj005'...
    };

% For each data set, there are 7 files
fnames = {...
    'epoch_conditions.mat' ...
    'ts_off_full.mat' ...
    'ts_off_left.mat' ...
    'ts_off_right.mat' ...
    'ts_on_full.mat' ...
    'ts_on_left.mat' ...
    'ts_on_right.mat' ...
    };

% Read / write the sample data
for s = 1:length(dirSubjects)
    
    thisdir = fullfile(savePth, dirSubjects{s});
    
    if ~exist(thisdir, 'dir'), mkdir(thisdir); end
    
    for f = 1:length(fnames)
        
        readPth  = fullfile(dirProject, dirSubjects{s}, fnames{f});
        
        writePth = fullfile(thisdir, fnames{f});
        
        fprintf('Downloading %s (please be patient).\n',fnames{f});
        
        urlwrite(readPth, writePth);
        
        fprintf('Downloading is done!\n');
        
    end
end

clear f files;
return



