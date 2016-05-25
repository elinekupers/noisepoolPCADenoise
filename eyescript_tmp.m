%% [Under construction] Script to get ms rate per condition, and plot mean rate over time.

% Note: Think about how putting these functions in our repository
toolbox_pth = '/Volumes/server/Projects/MEG/Eyetracking_scripts/';

% Add necessary paths:
addpath(fullfile(toolbox_pth));
addpath(genpath(fullfile(toolbox_pth,'toolboxes','mgl')));
addpath(genpath(fullfile(toolbox_pth,'toolboxes','mrToolsUtilities')));

% Get subject data
subjects = 6;
whichSubject = 6;
dataPath     = fullfile(dfdRootPath, 'exampleAnalysis', 'data');

% Get eye data
tmp          = load(sprintf(fullfile(dataPath, 'eye','s0%d_eyd.mat'),whichSubject));
eyd          = tmp.eyd; clear tmp;

% Get conditions
tmp          = load(sprintf(fullfile(dataPath, 's0%d_conditions.mat'),whichSubject));
conditions   = tmp.conditions;



startTime = NaN;
for ii = 1:length(eyd.messages)
    if strfind(eyd.messages(ii).message, 'MEG Trigger')
        startTime = ii; break;
    end
end


% Get eye data
% startTime   = esFindStart(eyd); % find the first message that contains MEG trigger
endTime     = 0; % number of messages to omit from end (if any)
timeLims    = [eyd.messages(startTime).time(1), eyd.messages(end-endTime).time];
blinkWinSec = [0.2 0.35];
s = msGetEyeData(eyd,timeLims,blinkWinSec);


% ------------- Delete irrelevant messages ----------------------------
eyd.messages = eyd.messages(1,startTime:end);

% ------------- Make a matrix for trigger nr and time stamp -----------
triggers = zeros(size(eyd.messages(1,:),2),2);

% ------------- Get trigger nr and time stamp -------------------------
for ii = 1:size(eyd.messages(1,:),2);                         % Last triggers while quiting the program are getting deleted
    triggers(ii,1) = str2num(eyd.messages(1,ii).message(14)); % Trigger nr
    triggers(ii,2) = eyd.messages(1,ii).time(1)-s.timeRaw(1); % Get timing and set time to zero
end

% ------------- Add triggers for the blank periods --------------------
onsets = ssmeg_trigger_2_onsets(triggers, whichSubject, 'eye');

% --------- Delete last 12 epochs since those are not recorded
% --------- % Question: Is this true for all subjects?
onsets = onsets(1:end-12);
conditions = conditions(1:end-12);

% Remove first epoch of flicker and non-flicker periods
% Define
badEpochs = zeros(size(onsets));
badEpochs(1:6:end) = 1;
% Remove
onsets = onsets(~badEpochs);
conditions = conditions(~badEpochs);


trialCount(whichSubject==subjects,:) = [size(find(conditions==3),1),size(find(conditions==1),1),size(find(conditions==5),1),size(find(conditions==7),1)];


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

% ------------- Define conditions -------------------------------------
blank   = sum(design,2)==0;
full	= design(:,1)==1;
left    = design(:,2)==1;
right	= design(:,3)==1;
conds   = {blank,full,left,right};
condsName = {'Blank','Full','Left','Right'};


vThres = 6;
msMinDur = 6;

% Concatenate all conditions for XY position and XY velocity
eyexyVel = [eyexVel(:), eyeyVel(:)];
eyexyPos = [eyexPos(:), eyeyPos(:)];



nn = 1;%:4
% Get velocity
thiseyexVel = eyexVel(:,conds{nn});
thiseyeyVel = eyeyVel(:,conds{nn});
% Get position
thiseyexPos = eyexPos(:,conds{nn});
thiseyeyPos = eyeyPos(:,conds{nn});

thiseyets = eyets(:,conds{nn});


% Concatenate all conditions for XY position and XY velocity
eyexyVel = cat(3,thiseyexVel, thiseyeyVel);
eyexyPos = cat(3,thiseyexPos, thiseyeyPos);

msVec = zeros(size(eyexyVel,2),1000);

for epoch = 1:size(eyexyVel,2);
    
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
%         s.xyPos(epoch)            = {squeeze(eyexyPos(:,epoch,:))};
        
        fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
            size(s.sacs{epoch},1), s.sacDetectRadius{epoch}(1), s.sacDetectRadius{epoch}(2));
        
        for msNr = 1:numSacs
            msVec(epoch,sacRaw(msNr,1)) = 1;
        end
    end
end

for epoch = 1:size(eyexyVel,2);
    msOnsets(epoch,:) = find(msVec(epoch,:));
end

