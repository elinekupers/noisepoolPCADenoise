function parsave(fname, varargin)
% Save variables to a file from within a parfor loop. A separate function
% is needed for this because the save function cannot be called directly
% from inside parfor loops
%
%  INPUTS:
%   fname: the full path of the matlab file to be saved
%   varagin: paired list {varname1, var1, varname2, var2, etc}
%
% Example:
%   fname = 'mymatfile.mat';
%   a = 1; b = 2;
%   parsave(fname, 'a', a, 'b', b);

if isodd(length(varargin))
    error('Varargin must contain paired entries, with a variable name followed by the variable')
end

% now we need to loop through varargin, and set each variable
for ii = 1:2:length(varargin)
    eval(sprintf('%s = varargin{%d};', varargin{ii}, ii+1))
end

% save
save(fname, varargin{1:2:end},'-v7.3');

end