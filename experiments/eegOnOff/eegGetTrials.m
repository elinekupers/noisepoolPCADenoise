function [pth, files, Axx, subj, thedate] = eegGetTrials(sessiondir, conditions)
% [pth, files, Axx] = eegGetTrials(sessiondir, conditions)


[subj, tmp] = strtok(sessiondir, '_');
thedate = strtok(tmp, '_');

root = strrep(which('setup.m'),'denoisesuite/setup.m','');
root = fullfile(root,'EEG','Data','PowerDiva');

d   = dir(sprintf('%s/%s/Exp*', root,sessiondir));
pth = sprintf('%s/%s/%s', root,sessiondir, d.name);

files = cell(1, length(conditions));
Axx   = cell(1, length(conditions));

for c = 1:length(conditions)
    
    % get filenames for condition c    
    files{c} = dir(fullfile(pth, sprintf('Raw_c%03d_*', conditions(c))));
    
        
    % load amplitude files for ON condition
    axxfile = fullfile(pth, sprintf('Axx_c%03d.mat', conditions(c)));
    if exist( axxfile, 'file'),  Axx{c} = load(axxfile);
    else                         Axx{c}.i1F1 = 18;    end
    
end
% 
% numfiles = cellfun(@numel, files);
% 
% % check that number of trials in all conditions is equal
% if ~isequal(min(numfiles), max(numfiles))
%     len = min(numfiles);
%     for c = 1:length(conditions)
%         files{c} = files{c}(end-len+1:end);        
%     end
% end
% 
% numfiles = cellfun(@numel, files);
% assert(isequal(min(numfiles), max(numfiles)))

return