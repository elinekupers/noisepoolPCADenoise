function savePth = DFDdownloaddata(savePth, whichDataSets, whichFiles)
% Download sample MEG data sets to be denoised by the 'Denoise Field Data'
% algorithm for the paper:
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% savePth = DFDdownloaddata(savePth, whichDataSets, whichFiles)
%
% Inputs
%   savePth: Path to store data. 
%                   [default = fullfile(DFDrootpath,'data')];
%   whichDataSets: Vector of data sets (1 to 8) 
%                   [default=1:8]
%   whichFiles: Vector of file numbers (there are 7 files per subject)
%                   [sefault = 1:7]
%
% Output
%   savePth: path where data was written
%
% Example 1: Download file 1 from subject 1
%   savePth = DFDdownloaddata([], 1, 1);
% Example 2: Download all files from subject 1
%   savePth = DFDdownloaddata([], 1);
% Example 3: Download all files from all subjects
%   savePth = DFDdownloaddata([]);

% Argument check
if notDefined('savePth'),       savePth       = fullfile(DFDrootpath, 'data'); end
if notDefined('whichDataSets'), whichDataSets = 1:8; end
if notDefined('whichFiles'),    whichFiles    = 1:7; end

% Site to retrieve the data
dirProject = 'http://psych.nyu.edu/winawerlab/denoiseFieldData';

% These are the subject data sets
dirSubjects = {...
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
for s = whichDataSets
    
    thisdir = fullfile(savePth, dirSubjects{s});
    
    if ~exist(thisdir, 'dir'), mkdir(thisdir); end
    
    fprintf('Downloading subject %d/%d (please be patient).\n',find(whichDataSets==s),length(whichDataSets));
    
    for f = whichFiles
        
        readPth  = fullfile(dirProject, dirSubjects{s}, fnames{f});
        
        writePth = fullfile(thisdir, fnames{f});
        
        urlwrite(readPth, writePth);
        
    end
end

fprintf('Downloading is done!\n');

return



