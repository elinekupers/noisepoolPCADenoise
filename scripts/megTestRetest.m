clear all;
inputDataDir = '/Volumes/HelenaBackup/denoisesuite/tmpmeg/';
outputFigDir = 'megfigs';
sensorDataStr = 'b2';
fitDataStr    = [sensorDataStr,'f_hpf2_fitfull30'];

sessionNums = [1,2,4,5];
maxperm = 100;
printFigsToFile = true;

%%
for k = 1:length(sessionNums)
    %% load fit data 
    sessionDir = megGetDataPaths(sessionNums(k));
    [resultsboot,r2boot] = megCombineSplithalvesHPC(fullfile(inputDataDir,'splitHalves'), sprintf('%s%s',sessionDir,fitDataStr),1:maxperm);
    % load original data 
    datafile = fullfile(inputDataDir,sprintf('%s%s',sessionDir,sensorDataStr));
    disp(datafile); load(datafile,'badChannels');
    
    %% plot data summary 
    npcs = resultsboot(1).opt.npcs;
    nchan = size(resultsboot(1).finalmodel.r2,2);
    
    xvaltrend = zeros(npcs+1,maxperm,2);
    optpcs    = zeros(maxperm,2);
    pcchan    = false(nchan,maxperm,2);
    noisepool = false(nchan,maxperm,2);
    beta      = zeros(3,nchan,maxperm,2);
    beta0     = zeros(3,nchan,maxperm,2);
    for np = 1:maxperm
        for nn = 1:2
            % top 10 channels
            pcchan(:,np,nn) = resultsboot(np,nn).pcchan{1};
            noisepool(:,np,nn) = resultsboot(np,nn).noisepool;
            % selected number of pcs
            optpcs(np,nn) = resultsboot(np,nn).pcnum(1);
            % curve for selecting PCs
            xvaltrend(:,np,nn) = mean(r2boot(:,pcchan(:,np,nn),nn,np),2);
            % get the beta values
            beta(:,:,np,nn) = resultsboot(np,nn).finalmodel.beta_md(:,:);
            beta0(:,:,np,nn) = resultsboot(np,nn).origmodel.beta_md(:,:);
        end
    end
    
    % prepare to plot 
    if k==1, fh1 = figure('Position',[1,500,1000,800]); end; figure(fh1);
    % Curve for selecting PCs, for first and second halves
    subplot(2,2,1); cla; hold on;
    plot(0:npcs,xvaltrend(:,:,1),'b');
    plot(0:npcs,xvaltrend(:,:,2),'r');
    xlabel('Number of PCs'); ylabel('R^2'); title('Top 10');
    makeprettyaxes;
    
    % Number of PCs chosen
    subplot(2,2,2); cla;
    hist(optpcs);
    optpcstat = prctile(optpcs,[16 50 84]);
    vline(optpcstat(2,1),'b');
    vline(optpcstat(2,2),'r');
    xlabel('Number of PCs'); ylabel('Frequency'); title('Chosen # PCs');
    makeprettyaxes;
    
    % Figure the channels that consistently show the biggest improvement
    subplot(2,2,[3,4]); cla; hold on;
    %plot(1:nchan,pcchan(:,:,1),'b');
    %plot(1:nchan,pcchan(:,:,2),'r');
    % top10 pcchannels
    pcchan2 = sum(reshape(pcchan,[nchan,maxperm*2]),2);
    plot(1:nchan,pcchan2,'k');
    % noisepool
    noise2 = sum(reshape(noisepool,[nchan,maxperm*2]),2);
    plot(1:nchan,noise2,'m');
    xlabel('Chan Number'); ylabel('# Perms'); xlim([1,nchan]); ylim([0,200]);
    legend('Top 10','Noise pool'); title('Channel selection frequency');
    makeprettyaxes;
    
    if printFigsToFile
        figurewrite(sprintf('sh%02d_%s%s', sessionNums(k), sessionDir, fitDataStr),[],[], sprintf('%s/s%d',outputFigDir,k), 1);
    else
        pause;
    end
    
    %% Print the beta estimation for the top 10 channels
    [sorted,sortedidx] = sort(pcchan2,'descend');
    chanNum = megGetOrigChannel(sortedidx(1:10)',badChannels,false);
    
    condNames  = {'FULL','RIGHT','LEFT'};
    if k== 1, fh2 = figure('Position',[1,500,400,800]); end; figure(fh2);
    
    for z = 1:10 % plot for each channel 
        % conds x maxperm x 2
        thisbeta  = squeeze(beta(:,sortedidx(z),:,:));
        thisbeta0 = squeeze(beta0(:,sortedidx(z),:,:));
        ctmp = [0.7,0.7,1; 1,0.7,0.7];
        
        for nn = 1:3 % for each condition 
            subplot(3,1,nn); cla; hold on;
            clear h hc h0 hc0
            for j = 1:2
                [h{j},hc{j}]=hist(thisbeta(nn,:,j),20);
                [h0{j},hc0{j}] = hist(thisbeta0(nn,:,j),20);
            end
            plot(hc{1},h{1},'b'); plot(hc{2},h{2},'r'); % after denoising
            plot(hc0{1},h0{1},'color',ctmp(1,:)); plot(hc0{2},h0{2},'color',ctmp(2,:)); % before denoising
            title(condNames{nn});
            makeprettyaxes;
            suptitle(sprintf('Chan %d (idx %d), Top10 freq=%0.1f',chanNum(z), sortedidx(z),sorted(z)/200));
        end
        if printFigsToFile
            figurewrite(sprintf('sh%02d_%s%s_beta%02d', sessionNums(k), sessionDir, fitDataStr,z),[],[], sprintf('%s/s%d',outputFigDir,k), 1);
        else
            pause;
        end
    end
end