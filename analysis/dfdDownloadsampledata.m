function savePth = dfdDownloadsampledata(savePth, whichSubjects, whichDataTypes)
% Download sample MEG data sets to be denoised by the 'Denoise Field Data'
% algorithm for the paper:
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% savePth = dfdDownloadSampleData(savePth, whichDataSets, whichFiles)
%
% Inputs
%   savePth: Path to store data. 
%                   [default = fullfile(dfdRootPath,'analysis','data')];
%
%   whichSubjects: Vector of one or more data sets (
%       [1:8]:          NYU datasets 1-8
%       [9:12] :        CiNet datasets 1-4
%       [13 15 17 19] : CiNet datasets 1-4 with SSS denoising
%       [14 16 18 20] : CiNet datasets 1-4 with tSSS denoising
%       [21:28] :       NYU datasets 1-8, with CALM denoising
%       [29:36] :       NYU datasets 1-8, with TSPCA denoising
%                   [default=1:8]
%
%   whichDataTypes: Cell array of one or more of  {'raw' ...
%                   'denoised 10 pcs' 'denoised 1-10 pcs' 'controls'}
%                   [default='raw']
%
% Output
%   savePth: path where data was written
%
% Example 1: Download raw data from subject 1
%   savePth = dfdDownloadSampleData([], 1, {'raw'});
% Example 2: Download raw data from all subjects
%   savePth = dfdDownloadSampleData();
% Example 1: Download denoised data from subject 1
%   savePth = dfdDownloadSampleData([], 1, {'denoised 10 pcs'});


% Argument check
if notDefined('savePth'),        savePth = fullfile(dfdRootPath, 'analysis', 'data'); end
if notDefined('whichSubjects'),  whichSubjects = 1:8; end
if notDefined('whichDataTypes'), whichDataTypes = {'raw'}; end

% Site to retrieve the data
dirProject = 'http://psych.nyu.edu/winawerlab/denoiseFieldData';

fnames = [];

for ii = 1:length(whichDataTypes)
    switch lower(whichDataTypes{ii})
        case 'raw'
            fnames = [fnames {'_conditions.mat' '_sensorData.mat'}];
        case 'denoised 10 pcs'
            fnames = [fnames {'_conditions.mat' '_denoisedData_bb.mat' '_denoisedData_sl.mat'}];
        case 'denoised 1-10 pcs'
            fnames = [fnames {'_conditions.mat' '_denoisedData_full_bb.mat' '_denoisedData_full_sl.mat'}];
        case 'controls'
            fnames = [fnames {'_denoisedData_control1_bb' '_denoisedData_control2_bb' '_denoisedData_control3_bb' '_denoisedData_control4_bb' '_denoisedData_control5_bb'}];
    end
end
fnames = unique(fnames);

% Read / write the sample data
for s = whichSubjects
        
    fprintf('Downloading subject %d .\n',s);
    
    for f = 1:length(fnames)
        
        fname = sprintf('s%02d%s', s, fnames{f});
        
        readPth  = fullfile(dirProject, fname);
        
        writePth = fullfile(savePth, fname);
        
        urlwrite(readPth, writePth);
        
    end
end

fprintf('Downloading is done!\n');

return



