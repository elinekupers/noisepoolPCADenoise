function savePth = dfdDownloadsampledata(savePth, whichDataSets)
% Download sample MEG data sets to be denoised by the 'Denoise Field Data'
% algorithm for the paper:
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% savePth = dfdDownloadsampledata(savePth, whichDataSets, whichFiles)
%
% Inputs
%   savePth: Path to store data. 
%                   [default = fullfile(DFDrootpath,'data')];
%   whichDataSets: Vector of data sets (1 to 8) 
%                   [default=1:8]
%
% Output
%   savePth: path where data was written
%
% Example 1: Download data from subject 1
%   savePth = dfdDownloadsampledata([], 1);
% Example 2: Download data from all subject
%   savePth = dfdDownloadsampledata();

% Argument check
if notDefined('savePth'),       savePth = fullfile(dfdRootPath, 'data', 'raw'); end
if notDefined('whichDataSets'), whichDataSets = 1:8; end

% Site to retrieve the data
dirProject = 'http://psych.nyu.edu/winawerlab/denoiseFieldData';


% For each data set, there are 7 files
fnames = {...
    '_conditions.mat' ...
    '_sensorData.mat' ...
    };

% Read / write the sample data
for s = whichDataSets
        
    fprintf('Downloading subject %d (please be patient).\n',s);
    
    for f = 1:length(fnames)
        
        fname = sprintf('s%02d%s', s, fnames{f});
        
        readPth  = fullfile(dirProject, fname);
        
        writePth = fullfile(savePth, fname);
        
        urlwrite(readPth, writePth);
        
    end
end

fprintf('Downloading is done!\n');

return



