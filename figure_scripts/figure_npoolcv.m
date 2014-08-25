%% THIS IS IN PROGRESS

inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/splitHalves';
fitDataStr   = 'b2fr_hpf2_fit10'; % fit data file string
sessionNums =  3;
maxperm = 100;

%%
for k = 1:length(sessionNums)
    % load fit data
    sessionDir = megGetDataPaths(sessionNums(k));
    filestr = fullfile(inputDataDir,sprintf('%s%s',sessionDir,fitDataStr));
    % combine across hpc results
    resultsboot = [];
    %idxboot1 = []; idxboot2 = [];
    for np = 1:maxperm
        currname = sprintf('%s_p%02d',filestr,np);
        matfile  = load(currname);
        resultsboot = cat(1,resultsboot,matfile.results);
        %idxboot1 = cat(2,idxboot1,[matfile.curr_idx1]);
        %idxboot2 = cat(2,idxboot2,[matfile.curr_idx2]);
    end
    for nt = 1:4
        % 3 conds x nchannels x 100 perms 
        origmodel = cat(1,resultsboot(:,nt).origmodel);
        origmodel_snr(:,:,:,nt) = cat(3,origmodel.beta_md)./cat(3,origmodel.beta_se);
        finalmodel = cat(1,resultsboot(:,nt).finalmodel);
        finalmodel_snr(:,:,:,nt) = cat(3,finalmodel.beta_md)./cat(3,finalmodel.beta_se);
    end
    
    % choose top 10
    finalsnr = mean(origmodel_snr(:,:,:,1),3);
    finalsnr = max(finalsnr);
    [~,idx] = sort(finalsnr,'descend');
    pcchan = false(size(resultsboot(1).noisepool));
    pcchan(idx(1:10))= 1;
    
    for icond = 1
        snrboot = squeeze(mean(origmodel_snr(icond,pcchan,:,:),2)); %nperms x 4
        snrbootint = prctile(snrboot,[16,50,84]);
        errorbar(1:4,snrbootint(2,:),snrbootint(1,:),snrbootint(3,:));
    end
    
end