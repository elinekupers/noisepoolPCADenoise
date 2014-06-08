function [resultsboot,r2boot] = megCombineSplithalvesHPC(resultsDir,filestr,perms)
% combines outputs of hpc script into single variables
% returns r2boot in the format of [npcs x nchannnels x nperm*2]
%         resultsboot in the format of nperm x 2 
% 
resultsboot = [];
r2boot      = [];
for np = perms
    currname = fullfile(resultsDir,sprintf('%s_p%02d',filestr,np));
    matfile  = load(currname);
    resultsboot = cat(1,resultsboot,matfile.results);
    r2boot      = cat(4,r2boot,matfile.r2);
end