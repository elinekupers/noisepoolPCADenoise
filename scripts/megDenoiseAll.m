function allResults = megDenoiseAll(sessionNums,dataDir,saveFlg,...
    dohpc,epochname,pcstop10,evalfunToCompute,sensorDataStr,saveDenoiseTs)
% denoise all sessions using commonly used parameters
% example: megDenoiseAll(1) : denoises session 1 and saves

% check input parameters
% ----------------------------------------------------------
if notDefined('dataDir'), dataDir = '/Volumes/server/Projects/MEG/GLMdenoised/tmpmeg'; end
if notDefined('saveFlg'), saveFlg = true; end
if notDefined('epochname') % whether to add an epoch group, if so, what the name is
    epochname = '';
    %epochname = '_epochGroup6o';
end
% whether to do high pass filtering as preprocessing
if notDefined('dohpc'),    dohpc = true; end
% whether to remove 10 pcs from the start
if notDefined('pcstop10'), pcstop10  = true; end
% which eval func to use ('bb','sl','bblog')
if notDefined('evalfunToCompute'), evalfunToCompute = {'bb'}; end
% which input sensor data set to use
if notDefined('sensorDataStr'), sensorDataStr = 'b2'; end
% save denoised time series
if notDefined('saveDenoiseTs'), saveDenoiseTs = false; end
% -----------------------------------------------------------

% load pre-saved frequencies
freq = load(fullfile(dataDir,'inputdata','megfreq'));
freq = freq.freq;
% create evoked function
evokedfun = @(x)getstimlocked(x,freq);
% create eval functions
slstr = '';
for ii = 1:length(evalfunToCompute)
    switch evalfunToCompute{ii}
        case 'bb'
            evalfun{ii} = @(x)getbroadband(x,freq);
            slstr = [slstr ''];
        case 'sl'
            evalfun{ii} = @(x)getstimlocked(x,freq);
            slstr = [slstr 'SL'];
        case 'bblog'
            evalfun{ii} = @(x)getbroadbandlog(x,freq);
            slstr = [slstr 'LG'];
    end
end

% -----------------------------------------------------------
% loop through sessions
for k = sessionNums
    % get and load input file
    [dataset,megDataDir] = megGetDataPaths(k,1:6);
    filename = [dataset,sensorDataStr];
    disp(filename)
    load(fullfile(dataDir,'inputdata',filename),'sensorData','design','badChannels');
    
    clear opt
    
    % figure out whether to add epoch groups
    % TODO here: can compute on the fly rather than use presaved ones
    if ~isempty(epochname)
        load(fullfile(dataDir,'epochGroups',[filename, epochname]));
    end
    % figure out whether to add preprocessing functions
    if dohpc
        opt.preprocessfun = @hpf;
        hpcstr = '_hpf2';
    else
        hpcstr = '';
    end
    
    % define options
    opt.freq = freq;   % frequencies
    opt.nboot = 100;  % how many times to bootstrap
    opt.resampling = {'boot','boot'}; % method for evokedfun glm and evalfun glm
    opt.npoolmethod = {'snr','n',75}; % noise pool selection method
    opt.pcselmethod = 'snr';          % pc selection method
    
    % if we decide to stop at 10 PCs
    if pcstop10
        opt.pcstop = -10;
        savestr = 'fit10';
    else
        savestr = 'fitful75';
    end
    
    % now denoise!
    if saveDenoiseTs
        [results,evalout] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    else
        [results,evalout,~,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    end
    
    % save stuff
    if saveFlg
        savename = fullfile(dataDir,sprintf('%sfr%s%s%s_%sp1k',filename,slstr,epochname,hpcstr,savestr));
        if saveDenoiseTs
            save(savename,'results','evalout','denoisedts');
        else
            save(savename,'results','evalout');
        end
        fprintf('data saved:%s\n', savename);
    end
    
    % return results
    allResults{k} = results;
end
