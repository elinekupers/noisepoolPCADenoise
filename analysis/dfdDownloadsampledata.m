function savePth = dfdDownloadsampledata(savePth, whichSubjects, whichDataTypes)
% Download sample MEG data sets to be denoised by the 'Denoise Field Data'
% algorithm for the paper:
%   Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)
%
% savePth = dfdDownloadsampledata(savePth, whichFiles)
% ----------------------------------------------------------------
% INPUTS:
% -----------------
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
%                   'denoised 10 pcs' 'denoised all pcs' 'controls'}
%                   [default='raw']
%
% OUTPUTS:
% -----------------
%   savePth: path where data was written
%
% EXAMPLES:
% -----------------
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
% dirProject = 'http://psych.nyu.edu/winawerlab/denoiseFieldData';
dirProject = 'https://osf.io';

urlStr = [];
fnames = [];

% Define the URL strings for the individual files
rawNYU =        {'ewtz9', 'gfhhk',  ...   conditions & data s01
                'jzx8u', '5yw3e',   ...   conditions & data s02
                'd8me2', 'zyrr2',   ...   conditions & data s03
                '92jhd', '5p8hh',   ...   conditions & data s04
                's7avt', 'hsbm2',   ...   conditions & data s05
                '7a5nq', 'nsqex',   ...   conditions & data s06
                'krdff', 'eycyu',   ...   conditions & data s07
                '6cbqr', 'tjj3v'};      % conditions & data s08


denoised10NYU = {'j6tzv', 'md8uu',  ...  denoised BB & SL s01
                'ffgev', 'fk8kb',   ...  denoised BB & SL s02
                'uu5b9', 'gw6up',   ...  denoised BB & SL s03
                'tex2u', '8aaxj',   ...  denoised BB & SL s04
                'ndw48', 'esvau',   ...  denoised BB & SL s05
                'xpc78', 'zh9sr',   ...  denoised BB & SL s06
                '8mrb5', 'sgeam',   ...  denoised BB & SL s07
                'uspmn', 'q8tpt'};     % denoised BB & SL s08

denoisedAllNYU = {'mf6nx','cw4fr',  ...   denoised BB & SL s01
                '3suq7', 'u77bm',   ...   denoised BB & SL s02
                '9hkns', '74vy5',   ...   denoised BB & SL s03
                '8m4vq', 'h53z8',   ...   denoised BB & SL s04
                's6m42', 'auzuh',   ...   denoised BB & SL s05
                'f4r6m', 'ach36',   ...   denoised BB & SL s06
                'hveh9', 'cpf9x',   ...   denoised BB & SL s07
                'gvpky', '9es6q'};      % denoised BB & SL s08

CALMNYU =       {'584hr', 'd92qf', ...    conditions & data s21
                'buqgg', '226gd', ...     conditions & data s22
                'g8cjm', 'xnmqq', ...     conditions & data s23
                '9y9ed', 'd92wr', ...     conditions & data s24
                '554gw', 'vnzjy', ...     conditions & data s25
                'kncfc', 'bqgz2', ...     conditions & data s26
                'yzegk', 'cp24s', ...     conditions & data s27
                'xtu4k', '8hgpy'};      % conditions & data s28

TSPCANYU =      {'6fwkd', 'wghxr',  ...   conditions & data s29
                '9jpyr', 'tj3cc',   ...   conditions & data s30
                'uw3ap', '2pexw',   ...   conditions & data s31
                'r5re8', 'qj4d2',   ...   conditions & data s32
                'th6xy', '8g3m4',   ...   conditions & data s33
                '8r97t', '6c947',   ...   conditions & data s34
                '6ufrt', 'y28mr',   ...   conditions & data s35
                'z5kq7', 'b3fs8'};      % conditions & data s36

controlNYU =    {'43p8e', 's7qam', 'suj6e', 'avqwf', 's2xsa', 'ugj6g',   ... Control 1-6 s01
                '6bfh6', '8xk5u', 'y4k2r', 'm77zu', '595cu', '364bq',    ... Control 1-6 s02
                'pevc7', '4r9xy','d827m', '5rw4y', 'kdfu8', 'gxgja',     ... Control 1-6 s03
                '3ctyg', '566es', '2f74r', '2wgs9', 'nfhea', 'ezr3e',    ... Control 1-6 s04
                'uz55w', 'accqa','msxjr','b39gy', 'n7nyg', 'zwgeb',      ... Control 1-6 s05
                'bj495', '4d5kr', 'cv4fp', 'vhtfy', '3nawg', 'wk38v',    ... Control 1-6 s06
                'pyvw5', 'swqhg', '4h8um', 'hfxz9', '8ykes', 'aj2hy',    ... Control 1-6 s07
                '9p9hy', 'h2czn', 'nn64m', '7udf5', '5jraa', 'wx3de'};     % Control 1-6 s08         

rawCiNET =      {'k3vub','n7j57',   ... % conditions & data s09
                'xdt9t', 'xz9d2',   ... % conditions & data s10
                'tj6wn', '9znqn',   ... % conditions & data s11
                'evkmk', 'suu3t'};      % conditions & data s12

TSSSCiNET =     {'7xj4j', '5vgh7',  ... conditions & data s14
                '', '',             ... left out SSS data (s15)
                'csk7e', 'u9xsy',   ... conditions & data s16
                '', '',             ... left out SSS data (s17)
                'w7d7n', 'fq4t2',   ... conditions & data s18
                '', '',             ... left out SSS data (s19)
                'fwvcj', 'ctc99'};      % conditions & data s20

% Concatenate the raw url's since we count subjects
raw = cat(2,rawNYU,rawCiNET,TSSSCiNET,CALMNYU,TSPCANYU);

% Get postfix of data file
for ii = 1:length(whichDataTypes)
    switch lower(whichDataTypes{ii})
        case 'raw'
            urlStr = [urlStr raw]; %#ok<*AGROW>
            fnames = [fnames, {'_conditions'},{'_sensorData'}];
        case 'denoised 10 pcs'
            urlStr = [urlStr denoised10NYU];
            fnames = [fnames, {'_denoisedData_bb'},{'_denoisedData_sl'}];
        case 'denoised all pcs'
            urlStr = [urlStr denoisedAllNYU];
            fnames = [fnames, {'_denoisedData_full_bb'},{'_denoisedData_full_sl'}];
        case 'controls'
            urlStr = [urlStr controlNYU];
            fnames = [fnames, {'_denoisedData_control'}];
    end
end
fnames = unique(fnames);

% Read / write the sample data
for s = whichSubjects
    
    fprintf('Downloading subject %d .\n',s);
    
    for f = 1:length(fnames)
        if fnames{f} ~= '_denoisedData_control'
            fname = sprintf('s%02d%s.mat', s, fnames{f});
            readPth  = fullfile(dirProject, urlStr{f}, '?action=download&version=1');
            writePth = fullfile(savePth, fname);
            websave(writePth,readPth);
        else
            for ncontrol = 1:6
                fname = sprintf('s%02d%s%d_bb.mat', s, fnames{f}, ncontrol);
                readPth  = fullfile(dirProject, urlStr{ncontrol+(6*(s-1))}, '?action=download&version=1');
                writePth = fullfile(savePth, fname);
                websave(writePth,readPth);
            end
        end
    end
end


fprintf('Downloading is done!\n');

return



