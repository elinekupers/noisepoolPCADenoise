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
dirPrefix = 'SSMEG_Dataset_';

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
    
    % Note that the first data set is called data set 4. The first 3 data
    % sets were pilots. For consistency with the naming scheme on our lab
    % server, we preserve the data set number.
    thisdir = sprintf('%s%02d', dirPrefix, s+3); 
    
    if ~exist(fullfile(savePth, thisdir), 'dir'), mkdir(fullfile(savePth, thisdir)); end
    
    fprintf('Downloading subject %d/%d (please be patient).\n',find(whichDataSets==s),length(whichDataSets));
    
    for f = whichFiles
        
        readPth  = fullfile(dirProject, thisdir, fnames{f});
        
        writePth = fullfile(savePth, thisdir, fnames{f});
        
        urlwrite(readPth, writePth);
        
    end
end

fprintf('Downloading is done!\n');

return



