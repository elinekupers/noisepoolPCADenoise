function dfdMakeFigure14()
%
% Function to get ms rate per condition, plot mean rate over time
% (directional histograms and bootstrapped distributions).

% dfdMakeFigure15()

% Eline Kupers, Helena X. Wang, Kaoru Amano, Kendrick N. Kay, David J.
% Heeger, Jonathan Winawer. (YEAR) Broadband spectral responses in visual
% cortex revealed by a new MEG denoising algorithm.
% (JOURNAL. VOLUME. ISSUE. DOI.)

% This script needs mgl and mrTools toolboxes!
% https://github.com/justingardner/mgl
% https://github.com/justingardner/mrTools
% 
% pth = '/Volumes/server/Projects/MEG/Eyetracking_scripts/toolboxes/mgl/'
% addpath(genpath(pth))

% pth = '/Volumes/server/Projects/MEG/Eyetracking_scripts/toolboxes/mrTools/'
% addpath(genpath(pth))


%% -------------------------------------
% ------------ Define Params -----------
% --------------------------------------

% ------------- Subject nr and folders -------------
whichSubject = 6; % Only subject 6,7,8 have eye tracking data
subjects     = whichSubject;
dataPath     = fullfile(dfdRootPath, 'analysis', 'data');
savePath     = fullfile(dfdRootPath, 'analysis', 'figures');
saveFigs     = false;
removeFirstEpoch = false;


% ------------- Experiment params -------------
condsName   = {'Blank','Full','Left','Right'};
condColors  = dfdGetColors(3);

% ------------- Eyetracking analysis params ---------
vThres = 6;      % Velocity threshold
msMinDur = 6;    % Microsaccade minimum duration (ms)

%% -------------------------------------
% ------- Load and prepare data --------
% --------------------------------------

% ------------- Get eye data -------------
tmp          = load(sprintf(fullfile(dataPath, 'eye','s0%d_eyd.mat'),whichSubject));
eyd          = tmp.eyd; clear tmp;

% ------------- Get conditions -------------
tmp          = load(sprintf(fullfile(dataPath, 's0%d_conditions.mat'),whichSubject));
conditions   = tmp.conditions;

% ------------- Get start time in eyelink data -------------
for ii = 1:length(eyd.messages)
    if strfind(eyd.messages(ii).message, 'MEG Trigger')
        startTime = ii; break;
    end
end

% ------------- Get messages -------------
endTime     = 0; % number of messages to omit from end (if any)
timeLims    = [eyd.messages(startTime).time(1), eyd.messages(end-endTime).time];
blinkWinSec = [0.2 0.35];
s           = msGetEyeData(eyd,timeLims,blinkWinSec);

% ------------- Delete irrelevant messages ----------------------------
eyd.messages = eyd.messages(1,startTime:end);

% ------------- Make a matrix for trigger nr and time stamp -----------
triggers = zeros(size(eyd.messages(1,:),2),2);

% ------------- Get trigger nr and time stamp -------------------------
for ii = 1:size(eyd.messages(1,:),2);                         % Last triggers while quiting the program are getting deleted
    triggers(ii,1) = str2num(eyd.messages(1,ii).message(14)); % Trigger nr
    triggers(ii,2) = eyd.messages(1,ii).time(1)-s.timeRaw(1); % Get timing and set time to zero
end

% ------------- Add triggers for the blank periods -------------
onsets = ssmeg_trigger_2_onsets(triggers, whichSubject, 'eye');

% ------------- Delete last 12 epochs since those are blank periods, which do not have triggers --------
% --------- % Question: Is this true for all subjects?
% onsets = onsets(1:end-12);
% conditions = conditions(1:end-12);

% ---------- Remove first epoch of flicker and non-flicker periods --------
badEpochs = zeros(size(onsets));
if removeFirstEpoch; badEpochs(1:6:end) = 1; end

% --------- Remove bad epochs ---------
onsets = onsets(~badEpochs);
conditions = conditions(~badEpochs);

% ------------- Check number of trials per conditions
trialCount(whichSubject==subjects,:) = [size(find(conditions==3),1),size(find(conditions==1),1),size(find(conditions==5),1),size(find(conditions==7),1)];


%% -------------------------------------
% ------------ Epoch data --------------
% --------------------------------------

% Convert x,y position and velocity vetors into epoched matrices (t x epoch)

% ------------- Define epochs in variables of eyetracking data --------
[eyets, ~]   = meg_make_epochs(s.time, onsets, [0 .999], 1000, 'eye');
[eyexPos, ~] = meg_make_epochs(s.xyPos(:,1), onsets, [0 .999], 1000, 'eye');
[eyeyPos, ~] = meg_make_epochs(s.xyPos(:,2), onsets, [0 .999], 1000, 'eye');
[eyexVel, ~] = meg_make_epochs(s.xyVel(:,1), onsets, [0 .999], 1000, 'eye');
[eyeyVel, ~] = meg_make_epochs(s.xyVel(:,2), onsets, [0 .999], 1000, 'eye');

% ------------- Make design matrix ------------------------------------
design = zeros(size(onsets,1),3);
design(conditions==1,1) = 1; % condition 1 is full field
design(conditions==5,2) = 1; % condition 5 is left field
design(conditions==7,3) = 1; % condition 7 is right field

% In case we dropped an epoch because we don't have enough timepoints to
% define the full last epoch
design = design(1:size(eyets,2),:);
conditions = conditions(1:size(eyets,2));

% ------------- Define conditions -------------------------------------
blank   = sum(design,2)==0;
full	= design(:,1)==1;
left    = design(:,2)==1;
right	= design(:,3)==1;
conds   = {blank,full,left,right};

% Concatenate all conditions for XY position and XY velocity
eyexyVel = [eyexVel(:), eyeyVel(:)];
eyexyPos = [eyexPos(:), eyeyPos(:)];

%% -------------------------------------
% ------- Eyetracking analysis ---------
% --------------------------------------


% For each condition
for nn = 1:numel(conds)
    
    % Get velocity
    thiseyexVel = eyexVel(:,conds{nn});
    thiseyeyVel = eyeyVel(:,conds{nn});
    % Get position
    thiseyexPos = eyexPos(:,conds{nn});
    thiseyeyPos = eyeyPos(:,conds{nn});
    % Get time
    thiseyets = eyets(:,conds{nn});
    
    % Concatenate all conditions for XY position and XY velocity
    eyexyVel = cat(3,thiseyexVel, thiseyeyVel);
    eyexyPos = cat(3,thiseyexPos, thiseyeyPos);
    
    numEpochs = size(eyexyVel,2);
    numTimePoints = size(thiseyexVel,1);
    msVec = zeros(numEpochs,numTimePoints);
    epochsWithMS = [];
    
    % For each second (==epoch)
    for epoch = 1:size(eyexyVel,2)
        
        % Define microsaccades
        [sacRaw,radius] = microsacc(squeeze(eyexyPos(:,epoch,:)),squeeze(eyexyVel(:,epoch,:)),vThres,msMinDur);
        
        % Remove the ones that occurr closely together (overshoot)
        numSacs = size(sacRaw,1);
        minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
        
        if numSacs ~= 0
            interSac = sacRaw(2:end,1)- sacRaw(1:end-1,2);
            sac = sacRaw([1; find(interSac > minInterSamples)+1],:);
            fprintf('%d rejected for close spacing\n', numSacs - size(sac,1));
            fprintf('%d saccades detected\n', size(sac,1));
            
            % Saved detected saccades into variable called 's'
            s.sacsRaw(epoch)          = {sacRaw};
            s.sacs(epoch)             = {sac};
            s.sacDetectRadius(epoch)  = {radius};
            s.eyeInfo.vThres(epoch,:)   = vThres;
            s.eyeInfo.msMinDur(epoch,:) = msMinDur;
            s.numSacs(epoch)            = numSacs;
            %         s.time(epoch)             = {thiseyets(:,epoch)};
            %             s.xyPos(epoch)            = {squeeze(eyexyPos(:,epoch,:))};
            
            fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
                size(s.sacs{epoch},1), s.sacDetectRadius{epoch}(1), s.sacDetectRadius{epoch}(2));
            
            
            
            
            
            
            % Make a binary microsaccade vector where each 1 is the start
            % of a ms in an epoch (1 s)
            for msNr = 1:numSacs
                msVec(epoch,sacRaw(msNr,1)) = 1;
            end
            
            epochsWithMS = [epochsWithMS epoch];
            
            % Get circular distribution of ms
            [rho{nn}{epoch},theta{nn}{epoch}] = msGetThetaRho(eyexyPos(:,epoch,:),s.sacs{epoch});
            
        end
    end
    
    % Concatenate ms vectors for each condition
    ms{nn} = msVec;
    
    epochNrInExperimentForThisCondition = find(conds{nn});
    allEpochsWithMS{nn} = epochNrInExperimentForThisCondition(epochsWithMS);
end

% Save vector with onset of microsaccades for each condition
save(sprintf(fullfile(dataPath, 'eye','s0%d_ms.mat'),whichSubject),'ms')
save(sprintf(fullfile(dataPath, 'eye','s0%d_epochswithms.mat'),whichSubject),'allEpochsWithMS')
%

%% -------------------------------------
% ------------ Get ms rates ------------
% --------------------------------------

% for epoch = 1:size(eyexyVel,2);
%     msOnsets(epoch,:) = find(msVec(epoch,:));
% end



% number of bins within one stimulus cycle for computing saccade rate
num_bins = 10;

% epoch time (seconds)
t = (1:numTimePoints)/numTimePoints;


% For simulation only
% make saccades by time-varying poisson process, where the rate is
% proportional to the phase of the stimulus
% r = ones(numEpochs,1) * (1+sin(2*pi*t*12))/100;
% r = r/mean(r(:))/1000;
% s = poissrnd(r);
% s = logical(s); % s indicates presence or absence of s
%

figure(1), clf, hold all; set(gcf, 'Color', 'w');
colors = [0 0 0; condColors];

for nn = 1:4
    
    thisMsVec = ms{nn};
    thisNrEpoch = size(thisMsVec,1);
    % stimulus phase (assumed to be 12 cycles per second)
    ph = mod(t*12*2*pi, 2*pi);
    ph = ones(thisNrEpoch,1)*ph;
    
    % get the stimulus phase for each saccade
    saccade_ph = ph(find(thisMsVec));
    
    % get the saccade count per time bin
    [y, bins] = hist(saccade_ph(:),num_bins);
    
    saccade_r = y / thisNrEpoch * num_bins;
    plot(bins/(2*pi)/12*1000,saccade_r, 'Color', colors(nn,:),'Marker','o', 'LineWidth',2)
end

legend(condsName,'Location', 'BestOutside');

ylabel('Saccade rate (per second)')
xlabel('Time (ms)')
title('Saccade rate in one stimulus period (one half cycle)')
makeprettyaxes(gca,18,18)

if saveFigs; hgexport(gcf, fullfile(savePath, sprintf('MSRateFlickerPeriod_subject%02d',whichSubject))); end

%% -------------------------------------
% --------- Bootstrap MS rates ---------
% --------------------------------------


figure(100); clf; set(gcf, 'Color', 'w'); hold all;
nboot  = 1000;
nbins  = [0:0.05:3];

condsName =  {'Blank','','Full','','Left','','Right'};

for nn = 1:4
    
    thisMsVec = ms{nn};
    epochnum = 1:size(thisMsVec,1);
    
    bootstat = bootstrp(nboot, @(x) mean(sum(thisMsVec(x, :),2)), epochnum);
    mdBootsStat = median(bootstat);
    y = hist(bootstat,nbins);
    plot(nbins, y, 'Color', colors(nn,:),'LineWidth',2);
    plot([mdBootsStat mdBootsStat],[0 500],'--','Color', colors(nn,:));
    
end

legend(condsName, 'Location', 'BestOutside')
xlabel('MS rate/s')
ylabel('Frequency')
makeprettyaxes(gca, 18,18);
title(sprintf('MS Distributions subject %02d', whichSubject));

if saveFigs; hgexport(gcf, fullfile(savePath, sprintf('MSDistributionEpoch_subject%02d',whichSubject))); end




% -------------------------------------
% ------ Circular distributions --------
% --------------------------------------

figure(101); clf; set(gcf, 'Color', 'w'); hold all;
radialLimit = [0 2];

for nn = 1:4;

    thisConditiontheta = theta{nn};
    thisConditionrho   = rho{nn};
    idx = cellfun(@isempty, thisConditiontheta);
    thisConditiontheta = thisConditiontheta(~idx);
    thisConditionrho = thisConditionrho(~idx);

    
    thetas = [];
    for ii = 1:numel(thisConditiontheta);
        thetas = vertcat(thetas, thisConditiontheta{ii});
    end
    
    rhos = [];
    for ii = 1:numel(thisConditionrho);
        rhos = vertcat(rhos, thisConditionrho{ii});
    end
   
    % plot individual microsaccade vectors
    subplot(2,4,nn);  
    h=polar2(thetas,rhos,radialLimit,'o'); 
    set(h,'markersize',2);
    set(h,'Color',colors(nn,:));
    set(h,'Color',colors(nn,:));
    title(sprintf('Microsaccade vectors %s',condsName{nn}))
    
    
    
    % circular distribution

    subplot(2,4,nn+4); cla;    
    t = 0:.01:(2*pi); 
    P1 = polar(t, 0.20 * ones(size(t))); hold all;
    P2 = polar(t, 0.15 * ones(size(t))); 
    P3 = polar(t, 0.10 * ones(size(t)));  
    P4 = polar(t, 0.05 * ones(size(t)));
%      set(P1, 'Visible', 'off')
    set(P1, 'Color', [.7 .7 .7])
    set(P2, 'Color', [.7 .7 .7])
    set(P3, 'Color', [.7 .7 .7])
    set(P4, 'Color', [.7 .7 .7])

    
    [tout, rout] = rose(thetas,(0:10:360)/180*pi);
    fh = polar(tout,rout/numel(theta{nn}));
    
%     rho_labels = {'0' '0.1' '0.2',};
%     rho_labels2 = {'0' '0.05' '.10', '.15' '0.2'};
%     
%     ff = findall(ax,'type','text');
%     t=strtrim(get(ff,'String'));
%     for r=1:length(rho_labels)
%     set(ff(strcmp(t,rho_labels{r})),'String',rho_labels2{r})
%     end
    
    set(fh, 'Color',colors(nn,:));
    set(fh, 'LineWidth', 2);
    title(sprintf('Angular distribution %s',condsName{nn}))
end


if saveFigs
    hgexport(gcf, fullfile(savePath, sprintf('MSVectorCircularDistributionEpoch_subject%02d',whichSubject))); 
end


