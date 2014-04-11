function [finalmodel,evalout,noisepool,denoisedspec] = denoisedata(design,data,evokedfun,evalfun,opt)

% first, get data dimensions 
[nchan,ntime,nepoch] = size(data); 
% handle inputs and options 
if notDefined('evokedfun'), evokedfun = @(x)getstimlocked(x,opt.freq); end
if notDefined('evalfun'),   evalfun   = @(x)getbroadband(x,opt.freq);  end
if notDefined('opt'),       opt       = struct(); end
if ~isfield(opt,'npoolmethod'), opt.npoolmethod = {'r2',[],'n',60};  end
if ~isfield(opt,'epochGroup'),  opt.epochGroup  = 1:nepoch;          end
if ~isfield(opt,'npcs'),        opt.npcs        = 30;                end
if ~isfield(opt,'verbose'),     opt.verbose     = true;              end
if ~isfield(opt,'xvalratio'),   opt.xvalratio   = -1;                end
if ~isfield(opt,'resampling'),  opt.resampling  = {'xval','xval'};   end
if ~isfield(opt,'pcstop'),      opt.pcstop      = 1.05;              end

if opt.verbose
    fprintf('---------------------------------------------------------------------\n');
    fprintf('(denoisedata) data dimenions: %d channels x %d time samples x %d epochs\n', nchan,ntime,nepoch);
    fprintf('---------------------------------------------------------------------\n');
end

% --------------------------------------------------------------
% perform fit to get R^2 values
% --------------------------------------------------------------
if opt.verbose, fprintf('(denoisedata) computing evoked model ...\n'); end
out = evalmodel(design,data,evokedfun,opt.resampling{1});

% --------------------------------------------------------------
% select noise pool
% --------------------------------------------------------------
if opt.verbose, fprintf('(denoisedata) selecting noise pool ...\n'); end
noisepool = selectnoisepool(out, opt.npoolmethod);
noisedata = data(noisepool,:,:);
if opt.verbose, fprintf('\t%d noise channels selected ...\n', sum(noisepool)); end
% check that the number of pcs we request isn't greater than the size of
% noisepool
if opt.npcs > sum(noisepool)
    fprintf('WARNING!!!: npcs > nnoise! Setting npcs to size of noise pool %d\n', sum(noisepool));
    opt.npcs = sum(noisepool);
end
    
% --------------------------------------------------------------
% compute PCs
% --------------------------------------------------------------
% pcs are stored in nrep cells; each cell is a matrix of [ntime x npcs]
nrep = max(opt.epochGroup);
if opt.verbose
    fprintf('(denoisedata) computing %d pcs for %d epoch groups ...\n', opt.npcs, nrep); 
end
pcs = cell(nrep,1);
for rp = 1:nrep
    % get current noise time series (ntime x nchan)
    currepochs = opt.epochGroup == rp;
    currnoise  = noisedata(:,:,currepochs);
    currnoise  = reshape(currnoise,[], ntime*sum(currepochs))';
    % unit-length normalize each time-series 
    temp = unitlengthfast(currnoise);
    % perform SVD and select top PCs
    %[u,s,v] = svd(temp);
    [coef,u,eigvals] = princomp(temp);
    u = u(:,1:opt.npcs);
    % scale so that std is 1 (ntime x npcs)
    pcs{rp} = bsxfun(@rdivide,u,std(u,[],1));
end

% --------------------------------------------------------------
% denoise in time and recompute spectral time series
% --------------------------------------------------------------
% evalaute each evalfun we pass in
if ~iscell(evalfun), evalfun = {evalfun}; end
nmodels = length(evalfun);
% loop through each pc
for p = 0:opt.npcs 
    if opt.verbose, fprintf('(denoisedata) denoising for %d pcs ...\n', p); end
    
    if p == 0
        denoiseddata = data;
    else
        denoiseddata = denoisetimeseries(data,pcs,p,opt.epochGroup);
    end
    
    % compute spectral time series and evaluate goodness of fit (by
    % applying evalfun we passed in)
    for fh = 1:nmodels
        [evalout(p+1,fh), denoisedst] = evalmodel(design,denoiseddata,evalfun{fh},opt.resampling{2},opt);
    end
end

% --------------------------------------------------------------
% % choose number of PCs, ala Kendrick
% --------------------------------------------------------------
% we probably eventually want to do something fancier, like fitting a curve
% Unlike fMRI data, EEG/MEG has fewer channels and curve is often not as
% smooth 
% QUESTION: how should we define xvaltrend?? in other words, how should R^2
% be combined across channels? now it's just an average across channels...
if opt.verbose, fprintf('(denoisedata) choosing pcs ...\n'); end
chosen = 0;
pcnum = zeros(1,nmodels);
for fh = 1:nmodels
    % npcs x channels, averaged across channels
    r2 = cat(1,evalout(:,fh).r2);
    xvaltrend = mean(r2,2);
    if opt.pcstop <= 0 % in this case, the user decides
        chosen = -opt.pcstop;
    else
        % this is the performance curve that starts at 0 (corresponding to 0 PCs)
        curve = xvaltrend - xvaltrend(1);
        % store the maximum of the curve
        mx = max(curve);
        % initialize (this will hold the best performance observed thus far)
        best = -Inf;
        for p=0:opt.npcs
            % if better than best so far
            if curve(1+p) > best
                % record this number of PCs as the best
                chosen = p;
                best = curve(1+p);
                % if we are within opt.pcstop of the max, then we stop.
                if best*opt.pcstop >= mx
                    break;
                end
            end
        end
    end
    % record the number of PCs
    pcnum(fh) = chosen;
    fprintf('\tevalfunc %d: %d pcs\n', fh, pcnum(fh));
end

% --------------------------------------------------------------
% get final model and denoised time series
% --------------------------------------------------------------
% pull out the final models and add pcnum to it
for fh = 1:nmodels
    denoiseddata = denoisetimeseries(data,pcs,pcnum(fh),opt.epochGroup); 
    [finalmodel(fh), denoisedspec{fh}] = evalmodel(design,denoiseddata,evalfun{fh},'full',opt);
end
for fh = 1:nmodels, finalmodel(fh).pcnum = pcnum(fh);end

% % save denoised time series, if requested
% if nargout>3
%     denoisedspec = cell(1,nmodels);
%     for fh = 1:nmodels
%         denoisedspec{fh} = denoisetimeseries(data,pcs,pcnum(fh),opt.epochGroup); 
%     end
% end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out,datast] = evalmodel(design,data,func,how,opt)
% INPUTS:
% data   : [channels x time x epochs]
% design : [epochs x nvectors]
% func   : a function that summarizes data into [epochs x channels]
% how    : defines how fitting is done ['full' or 'xval'(default)]
%
% OUTPUTS:
% r2     : goodness of fit   [channels x 1]
% beta   : betas for the fit [nperms x channels]

% check inputs
if notDefined('how'),  how  = 'xval';   end
if notDefined('opt'),  opt  = struct(); end
if ~isfield(opt,'verbose'),    opt.verbose    = true; end
if ~isfield(opt,'maxpolydeg'), opt.maxpolydeg = 0;     end
if ~isfield(opt,'xvalratio'),  opt.xvalratio  = -1;    end
if ~isfield(opt,'xvalmaxperm'),opt.xvalmaxperm =500;   end

% datast should be dimensions [epochs x channels]
datast = func(data);
nepochs = size(datast,1);
% sanity checks
assert(nepochs==size(design,1));
assert(sum(isnan(datast(:)))==0 && sum(isinf(datast(:)))==0);

% remove polynomial from data and design
% for now this is usually just a constant term, so we're just de-meaning
pmatrix = constructpolynomialmatrix(nepochs,0:opt.maxpolydeg);
pmatrix = projectionmatrix(pmatrix);
datast  = pmatrix*datast;
design  = pmatrix*design;

switch how
    case 'full'
        if opt.verbose, fprintf('\tfull fit: no resampling\n'); end
        % do glm
        beta = design \ datast;
        modelfit = design*beta;
        % compute goodness of fit
        r2 = calccod(modelfit,datast,1,[],0);
        
        % save into output struct
        out = struct('r2',r2,'beta',beta,'modelfit',modelfit);
        
    case 'xval'
        % figure out how to choose train and test data
        if opt.xvalratio == -1 % do n-fold cross validation
            epochs_test = (1:nepochs)';
            if opt.verbose, fprintf('\txval: n-fold/leave-on-out\n'); end
        else % else, divide up according to opt.xvalratio
            ntest = round(nepochs * opt.xvalratio);
            ntest(ntest<1) = 1; ntest(ntest>=nepochs)=nepochs-1;
            % if too many permutations, so we choose randomly
            % just an arbitrary definition for "too many"
            if nchoosek(nepochs,ntest) > max(nepochs*3, opt.xvalmaxperm)
                epochs_test = zeros(opt.xvalmaxperm,ntest);
                for nn = 1:opt.xvalmaxperm
                    tmp = randperm(nepochs);
                    epochs_test(nn,:) = tmp(1:ntest);
                end
            else
                epochs_test = nchoosek(1:nepochs,ntest);
            end
            if opt.verbose
                fprintf('\txval: ntest %d; ntrain %d (%0.2f test)\n', ntest, nepochs-ntest,opt.xvalratio); 
            end
        end
        
        % for each train/test permutation
        modelfit = []; beta = []; r2perm = [];
        for nn = 1:size(epochs_test,1)
            curr_test = epochs_test(nn,:);
            curr_train= setdiff(1:nepochs,curr_test);
            % do glm on training data
            % beta_train = [n x channels]
            beta_train= design(curr_train,:)\ datast(curr_train,:);
            modelfit_test = design(curr_test,:)*beta_train;
            % save prediction for this epcoh
            modelfit = cat(1,modelfit,modelfit_test); %[epochs x channels]
            beta     = cat(3,beta,    beta_train);    %[n x channels x perms]
            r2perm   = cat(1,r2perm,  calccod(modelfit_test,datast(curr_test,:),[],0)); % [perms x 1]
        end
        % both data and predicted are [epochs x channels]; r2 = [1 x channels]
        r2 = calccod(modelfit,datast(vectify(epochs_test'),:),[],0);
        
        % save into output struct
        out = struct('r2',r2,'beta',beta,'modelfit',modelfit, 'epochs_test', epochs_test, 'r2perm', r2perm);
    
    otherwise
        error('(denoisedata:evalmodel) resampling method not parsed: %s', how);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function noisepool = selectnoisepool(out,npoolmethod)
% selects noise channels
% npoolmethod is a cell array
%   entry 1 defines what variable is used ['r2'(default) | 'beta']
%   entry 2 (optional) defines function on out [default [] ] 
%           can pass in a function to calculate snr, for example
%   entry 3 defines how noise channels are selected ['n' (default)|'thres']
%           if 'n', then we choose the worst X channels 
%           if 'thres', then we choose channels with entry1 < X
%   entry 4 defines X 

% check inputs 
var  = npoolmethod{1}; if isempty(var), var = 'r2'; end
if ~strcmp(var,'r2') && ~strcmp(var,'beta')
    error('(denoisedata:selectnoisepool): npoolmethod{1} not parsed: %s', var);
end
func = npoolmethod{2};
compmethod = npoolmethod{3}; if isempty(compmethod), compmethod = 'n'; end
compval    = npoolmethod{4}; 
if isempty(compval) 
    if strcmp(compmethod,'n'), compval = 60;
    elseif strcmp(compmethod,'thres') compval = 0; end
end

% compute 
vals = out.(var);
if ~isempty(func), vals = func(vals); end
if strcmp(compmethod,'n')
    noisepool = false(1,size(out.r2,2));
    [~,sortedinds] = sort(vals,'ascend');
    noisepool(sortedinds(1:compval)) = true;
end
if strcmp(compmethod,'thres')
    noisepool = vals < compval;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function denoiseddata = denoisetimeseries(data,pcs,p,epochGroup)
% denoise time series
% INPUTS:
% data   : [channels x time x epochs]
% pcs    : nrep cells; each is a matrix of [ntime x npcs]
% p      : how many number of pcs we want to use 
% how    : defines how epochs are grouped 
%
% OUTPUTS:
% denoiseddata : [channels x time x epochs]
%
[nchan,ntime,nepoch] = size(data);
nrep = max(epochGroup);

denoiseddata = zeros(nchan,ntime,nepoch);
cummnepoch   = 0;
% project out appropriate number of pcs from each epoch or epoch group
for rp = 1:nrep
    % get time series (ntime x nchan) for current epoch group
    currepochs = epochGroup == rp;
    currnepoch = sum(currepochs);
    currsig    = data(:,:,currepochs);
    currsig    = reshape(currsig,[nchan,ntime*currnepoch])';
    
    % denoise data for this epoch and this number of pcs
    currdenoisedsig = currsig - pcs{rp}(:,1:p)*(pcs{rp}(:,1:p)\currsig);
    
    % reshape into ntime x nepoch x nchan
    currdenoisedsig = reshape(currdenoisedsig,[ntime, currnepoch, nchan]);
    % save into denoiseddata for this number of pcs [nchan x ntime x nrep]
    denoiseddata(:,:,cummnepoch+(1:currnepoch)) = permute(currdenoisedsig, [3,1,2]);
    cummnepoch = cummnepoch + currnepoch;
    
    % sanity check
    if rp==nrep, assert(cummnepoch(end) == nepoch); end
end
