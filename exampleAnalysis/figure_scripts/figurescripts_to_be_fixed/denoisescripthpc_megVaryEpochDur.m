function denoisescripthpc_megVaryEpochDur()
% denoise as a function of denoising epoch duration 
% you can modify the script take parameter k for session number

root = '/scratch/hxw201/denoisedir/';
inputdir = 'inputdata';
addpath(fullfile(root,'aux'));
addpath(fullfile(root,'funcs'));
datapaths = sessions();

dohpc     = true;

freq = load(fullfile(root,inputdir,'megfreq'));
freq = freq.freq;
%warning off;
%epochDurs = [1,3*2.^(0:7),inf];
epochDurs = [1,3,6,12,24,36,72,inf];
npcs      = [5,10:10:70];

for k = [3:6,9:12]%10:12%[3:6,9:12]
    filename = [datapaths{k},'b2'];
    disp(filename);
    %filename = datapaths{k};
    load(fullfile(root,inputdir,filename),'sensorData','design','badChannels');
    load(fullfile(root,inputdir,'okEpochs',filename),'okEpochs');
    if sum(okEpochs) ~= size(sensorData,3)
        error('wrong okEpochs');
    end
    
    clear opt
    
    if dohpc
        opt.preprocessfun = @hpf;
        hpcstr = '_hpf2';
    else
        hpcstr = '';
    end
    
    evokedfun = @(x)getstimlocked(x,freq);
    %evalfun   = {@(x)getbroadband(x,freq), @(x)getstimlocked(x,freq)};
    %evalfun = {@(x)getbroadbandlog(x,freq)};
    evalfun = {@(x)getbroadband(x,freq), @(x)getbroadbandlog(x,freq)};
    for kk = 1:length(evalfun), evalfunstr{kk} = func2str(evalfun{kk}); end
    opt.freq = freq;
    %opt.npcs2try    = 70;
    %opt.xvalmaxperm = 100;
    opt.resampling  = {'boot','boot'};
    opt.nboot = 100;
    %opt.resampling = {'xval','xval'};
    %opt.fitbaseline = true;
    %opt.savepcs = true;
    opt.pcselmethod = 'snr';
    
    opt.npoolmethod = {'snr','n',75};
    %opt.npoolmethod = {'r2','n',75};
    %opt.npoolmethod = {'r2','thres',0};
    
    allResults = [];
    for ii = 1:length(epochDurs) % iterate through epoch duration fo doing PCA
        
        if isinf(epochDurs(ii))
            epochGroup = ones(size(sensorData,1),1);
        else
            epochGroup = megEpochGroup(okEpochs,epochDurs(ii),0);
        end
        opt.epochGroup = epochGroup;
        clear results;
        fprintf('epochDur = %d; %d epochs\n', epochDurs(ii), max(epochGroup));
        
        for jj = 1:length(npcs) % iterate through number of PCs projected out 
            opt.pcstop = -npcs(jj);
            fprintf('\tnpcs = %d\n', npcs(jj));
            
            if jj == 1
                [results] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
                noisepooldef = results.noisepool;
            else
                [results] = denoisedata(design,sensorData,noisepooldef,evalfun,opt);
            end
            allResults{ii,jj} = results;
        end
    end
    
    %[results,evalout] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    savename = fullfile(root,'data',sprintf('%sfr%s_fitfull75p1k_varyEpochs',filename,hpcstr));
    save(savename,'allResults','evalfunstr','epochDurs','npcs');
    %save(savename,'results','evalout','denoisedts');
    fprintf('data saved:%s\n', savename);
    clear okEpochs
end
%warning on;