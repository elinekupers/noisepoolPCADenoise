tmpmegdir   = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
sessionNums = 3;
k = 1;

sessionDir = megGetDataPaths(sessionNums(k));
ppfit = 'b2f_hpf2_fitfull30';
maxperm = 100;
[resultsboot,r2boot] = megCombineSplithalvesHPC(fullfile(tmpmegdir,'splitHalves'), sprintf('%s%s',sessionDir,ppfit),1:maxperm);

%%
npcs = resultsboot(1).opt.npcs;
xvaltrend = zeros(npcs+1,maxperm,2);
optpcs    = zeros(maxperm,2);
beta      = zeros(3,size(resultsboot(1).finalmodel.r2,2),maxperm,2);
for np = 1:maxperm
    for k = 1:2
        pcchan = resultsboot(np,k).pcchan{1};
        optpcs(np,k) = resultsboot(np,k).pcnum(1);
        xvaltrend(:,np,k) = mean(r2boot(:,pcchan,k,np),2);
        beta(:,:,np,k) = resultsboot(np,k).finalmodel.beta_md(:,:);
    end
end
hold on;
plot(0:npcs,xvaltrend(:,:,1),'b');
plot(0:npcs,xvaltrend(:,:,2),'r');

%%

