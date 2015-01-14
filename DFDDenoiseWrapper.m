function allResults = DFDDenoiseWrapper(sessionNums,dataDir,saveFlg,...
    dohpc,epochname,pcstop10,evalfunToCompute,sensorDataStr,saveDenoiseTs)
%% Prepare all requested sessions and denoise with predefined parameters
%  by the 'Denoise Field Data' algorithm for the paper:
%
%   AUTHORS. YEAR. TITLE. JOURNAL. VOLUME. ISSUE. DOI.
%
% allResults = DFDDenoiseWrapper(sessionNums,dataDir,saveFlg, ...
%    dohpc,epochname,pcstop10,evalfunToCompute,sensorDataStr,saveDenoiseTs)
%
% Inputs
%   sessionNums:        Vector of data sets (1 to 8)
%                       [default=1:8]
%   dataDir:            String to define data directory
%                       [default=ullfile(DFDrootpath, 'data')]
%   saveFlg:            Boolean whether to save preprocessed datasets or not
%                       [default=true]
%   dohpc:              Boolean whether to highpass filter or not
%                       [default=true]
%   epochname:          String whether to add an epoch group, if so, what the name is
%                       [default='']
%   pcstop10:           Boolean whether to jump to 10 pcs from the start
%                       and don't evaluate PCs in between.
%                       [default=true]
%   evalfunToCompute:   cell to define which function to use to evaluate 
%                       [default={'bb'}]
%   sensorDataStr:      String to define which input sensor data set to use
%                       [default='b2']
%   saveDenoiseTs:      Boolean whether to save denoised time series or not
%                       [default=false]
%
% Outputs
%   allResults:         A cell containg several results:
%                           - 'results' (...)
%                           - 'evalout' (output of specified evaluation function) 
%                           - 'denoisedts' (if requested)
%
% Example: Denoise broadband of session 1 and save data
%   allResults = DFDDenoiseWrapper(1)
% 
% 


%% check input parameters
% ----------------------------------------------------------
if notDefined('dataDir'), dataDir = fullfile(DFDrootpath, 'data'); end
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
%% -----------------------------------------------------------

% load pre-saved frequencies, we use a 1 second epoch and go up to 150 Hz
% and have a 12 Hz stimulus locked signal.
T = 1; fmax = 150; slF = 12;
freq = megGetSLandABfrequencies((0:fmax)/T, T, slF/T);


% ---------- I replaced these lines with the two above ----------------
% freq = load('megfreq'); %% <--- HACK: I just put the file in the folder, but we need to find where Helena defined this
% freq = freq.freq;

%% ----------------------------------------------- 

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

%% -----------------------------------------------------------
% loop through sessions
for k = sessionNums
    % get and load input file
    [dataset,megDataDir] = DFDgetdatapaths(k,1:6,dataDir);
    filename = [dataset,sensorDataStr];
    disp(filename)
    load(fullfile(dataDir,'savedProcData',filename),'sensorData','design','badChannels');
    
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
    opt.nboot = 1000;  % how many times to bootstrap
    opt.resampling = {'boot','boot'}; % method for evokedfun glm and evalfun glm
    opt.npoolmethod = {'snr','n',75}; % noise pool selection method
    opt.pcselmethod = 'snr';          % pc selection method
    
    % if we decide to stop at 10 PCs
    if pcstop10
        opt.pcstop = -10;
        savestr = 'fit10';
    else
        savestr = 'fitfull75';
    end
    
    % now denoise!
    if saveDenoiseTs
        [results,evalout,~,denoisedts] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    else
        [results,evalout] = denoisedata(design,sensorData,evokedfun,evalfun,opt);
    end
    
    % save stuff
    if saveFlg
        savename = fullfile(dataDir,'savedProcData',sprintf('%sfr%s%s%s_%sp1k',filename,slstr,epochname,hpcstr,savestr));
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
