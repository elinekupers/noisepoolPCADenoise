function dfdEyeScript(subjects)

%
% dfdEyeScript(subjects)
%
% INPUTS:
% subjects    : Number of subjects one would like to denoise
%
% DESCRIPTION: Function to analyze eyetracking data recorded during MEG visual steady
% state experiment, for subject 6, 7 and 8.
%
% DEPENDENCIES: This function depends on function from the meg_utils
% repository
%
% AUTHORS. YEAR. TITLE. JOURNAL.

% =========================================================================
% =============== Check options and load in toolboxes =====================
% =========================================================================

% Check input
if subjects < 6; disp('Subject does not have eye tracking data'); end;

% Note: Think about how whether putting these functions in our repository
toolbox_pth = '/Volumes/server/Projects/MEG/Eyetracking_scripts/';


% Add necessary paths:
addpath(fullfile(toolbox_pth));
addpath(genpath(fullfile(toolbox_pth,'toolboxes','mgl')));
addpath(genpath(fullfile(toolbox_pth,'toolboxes','mrToolsUtilities')));
addpath(genpath('~/Applications/mgl/mgllib'));

% Check options:
saveEyd           = true;  % Convert edf to eyd.mat file and save it?
saveFigures       = true;  % Save images?
deleteFirstLast   = false; % Delete first and last epoch?
saveData          = true;  % Save data?

% =========================================================================
% ================= Define paths and load in data =========================
% =========================================================================
for whichSubject = subjects
    edffile = dir(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye','s0%d_eyelink.edf'),whichSubject));
    tmp = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 's0%d_conditions.mat'),whichSubject)); conditions = tmp.conditions;
    
    
    if saveEyd
        eyd = mglEyelinkEDFRead(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye',edffile.name)));
        save(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye', 's0%d_eyd.mat'),whichSubject));
    else
        eyd = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye','s0*_eyd.mat'),whichSubject));
    end
    
    % =====================================================================
    % ================ Get eye traces (blinks are removed) ================
    % =====================================================================
    
    % Figure out the start and end time of your experiment, which will be a
    % an argument for the msGetEyeData function.
    % Blink window is pretty conservative right now. could be [-0.1,0.1] (100
    % ms either way)
    startTime   = 8;
    endTime     = 0;
    timeLims    = [eyd.messages(startTime).time(1), eyd.messages(end-endTime).time];
    blinkWinSec = [0.2 0.35];
    s = msGetEyeData(eyd,timeLims,blinkWinSec);
    
    % =====================================================================
    % ===================== Get triggers and timings ======================
    % =====================================================================
    
    % Use MEG messages to get the eyemovements for separate conditions
    
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
    
    % ------------- Define epochs in variables of eyetracking data --------
    [eyets, ~]   = meg_make_epochs(s.time, onsets, [0 .999], 1000);
    [eyexPos, ~] = meg_make_epochs(s.xyPos(:,1), onsets, [0 .999], 1000);
    [eyeyPos, ~] = meg_make_epochs(s.xyPos(:,2), onsets, [0 .999], 1000);
    [eyexVel, ~] = meg_make_epochs(s.xyVel(:,1), onsets, [0 .999], 1000);
    [eyeyVel, ~] = meg_make_epochs(s.xyVel(:,2), onsets, [0 .999], 1000);
        
    % ------------- Make design matrix ------------------------------------
    design = zeros(size(eyets,2),3);
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
    
    
 
    %% ====================================================================
    %  ============ Plot eye traces for visual inspection =================
    %  ====================================================================

    %% Plot X COORDINATES
    deglims = 20;
    figure(1); clf; set(gcf,'Color', 'w');
    subplot(2,2,1);
    
    % All eyetracking data
    plot(s.time/1000,s.xyPos(:,1), 'k'); hold on;
     % Plot per stimulus condition
    for nn = 1:4
         plot(eyets(:,conds{nn}),eyexPos(:,conds{nn}));
    end
    grid on;
    ylim(deglims*[-1,1]); xlabel('Time (s)'); ylabel('X (deg)');
    
    
    %% Plot Y COORDINATES  
    subplot(2,2,3);
    
    % All eyetracking data
    plot(s.time/1000,s.xyPos(:,2),'k'); hold on;    
    % Plot per stimulus condition
    for nn = 1:4
         plot(eyets(:,conds{nn}),eyeyPos(:,conds{nn}));
    end
    
    grid on;
    ylim(deglims*[-1,1]); xlabel('Time (s)'); ylabel('Y (deg)');
    
    %% Plot XY ON GRID
    subplot(2,2,[2,4]);

    % All eyetracking data
    plot(s.xyPos(:,1),s.xyPos(:,2), 'k'); axis square; grid on; hold on;
    % Plot per stimulus condition
    for nn = 1:4
         plot(eyexPos(:,conds{nn}),eyeyPos(:,conds{nn}));
    end
     
    xlim(deglims*[-1,1]); ylim(deglims*[-1,1]);
    xlabel('X (deg)', 'Fontsize', 20);  ylabel('Y (deg)', 'Fontsize', 20);
    set(gca, 'FontSize', 20);

    
    %% =================================================
    %  =============== Detect saccades =================
    %  =================================================
    
    % Eyelink also detects saccades, but the algorithm isn't that reliable.
    % This is only if you actually want to do analyses with the saccades (e.g.,
    % look at how microsaccades are modulated by attention). If just for
    % monitoring goodness of fixation, you can probably just rely on the
    % Eyelink detected saccades (eyd.saccades)
    vThres = 6;
    msMinDur = 6;
    
    % Concatenate all conditions for XY position and XY velocity
    eyexyVel = [eyexVel(:), eyeyVel(:)];
    eyexyPos = [eyexPos(:), eyeyPos(:)];
    
    %% ALL data
    [sacRaw,radius] = microsacc(s.xyPos,s.xyVel,vThres,msMinDur);
    
    % remove the ones that occurr closely together (overshoot)
    numSacs = size(sacRaw,1);
    minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
    interSac = sacRaw(2:end,1)- sacRaw(1:end-1,2);
    sac = sacRaw([1; find(interSac > minInterSamples)+1],:);
    fprintf('%d rejected for close spacing\n', numSacs - size(sac,1));
    
    fprintf('%d saccades detected\n', size(sac,1));
    
    % saved detected saccades into s
    s.sacsRaw          = sacRaw;
    s.sacs             = sac;
    s.sacDetectRadius  = radius;
    s.eyeInfo.vThres   = vThres;
    s.eyeInfo.msMinDur = msMinDur;
    
    fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
        size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));
    
    % look at detected saccades
    figure(2); clf; set(gcf,'Color', 'w');
    msSacStats1(s);
    title('Angular distribution ALL');
    
    %% Epoched data
    
    [sacRaw,radius] = microsacc(eyexyPos,eyexyVel,vThres,msMinDur);
    
    % remove the ones that occurr closely together (overshoot)
    numSacs = size(sacRaw,1);
    minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
    
    %% DEBUG FROM HERE %%
    interSac = sacRaw(2:end,1)- sacRaw(1:end-1,2);
    sac = sacRaw([1; find(interSac > minInterSamples)+1],:);
    fprintf('%d rejected for close spacing\n', numSacs - size(sac,1));
    
    fprintf('%d saccades detected\n', size(sac,1));
    
    for nn = 1:4
        % saved detected saccades into s
        s.sacsRaw          = sacRaw(:,conds{nn});
        s.sacs             = sac(:,conds{nn});
        s.sacDetectRadius  = radius(:,conds{nn});
        s.eyeInfo.vThres   = vThres;
        s.eyeInfo.msMinDur = msMinDur;
        s.time             = eyets(:,conds{nn});
        s.xyPos            = eyexyPos(:,conds{nn});

        fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
            size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));

        % look at detected saccades
        figure(2+nn); clf; set(gcf,'Color', 'w');
        msSacStats1(s);
        title(sprintf('Angular distribution %s',condsName{nn}));
    end
    
    %% OFF
    
    [off_sacRaw,off_radius] = microsacc(off_xypos,off_xyvel,vThres,msMinDur);
    
    % remove the ones that occurr closely together (overshoot)
    off_numSacs = size(off_sacRaw,1);
    minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
    off_interSac = off_sacRaw(2:end,1)- off_sacRaw(1:end-1,2);
    off_sac = off_sacRaw([1; find(off_interSac > minInterSamples)+1],:);
    fprintf('%d rejected for close spacing\n', off_numSacs - size(off_sac,1));
    
    fprintf('%d saccades detected\n', size(off_sac,1));
    
    % saved detected saccades into s
    s.sacsRaw = off_sacRaw;
    s.sacs = off_sac;
    s.sacDetectRadius  = off_radius;
    s.eyeInfo.vThres   = vThres;
    s.eyeInfo.msMinDur = msMinDur;
    s.time = eye_ts_off_epoched;
    s.xyPos = off_xypos;
    
    fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
        size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));
    
    % look at detected saccades
    figure(4); clf; set(gcf,'Color', 'w');
    msSacStats1(s);
    title('Angular distribution OFF');
    
    %% LEFT
    
    [left_sacRaw,left_radius] = microsacc(left_xypos,left_xyvel,vThres,msMinDur);
    
    % remove the ones that occurr closely together (overshoot)
    left_numSacs = size(left_sacRaw,1);
    minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
    left_interSac = left_sacRaw(2:end,1)- left_sacRaw(1:end-1,2);
    left_sac = left_sacRaw([1; find(left_interSac > minInterSamples)+1],:);
    fprintf('%d rejected for close spacing\n', left_numSacs - size(left_sac,1));
    
    fprintf('%d saccades detected\n', size(left_sac,1));
    
    % saved detected saccades into s
    s.sacsRaw = left_sacRaw;
    s.sacs = left_sac;
    s.sacDetectRadius  = left_radius;
    s.eyeInfo.vThres   = vThres;
    s.eyeInfo.msMinDur = msMinDur;
    s.time = eye_ts_on_left_epoched;
    s.xyPos = left_xypos;
    
    fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
        size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));
    
    % look at detected saccades
    figure(5); clf; set(gcf,'Color', 'w');
    msSacStats1(s);
    title('Angular distribution LEFT');
    
    %% RIGHT
    
    [right_sacRaw,right_radius] = microsacc(right_xypos,right_xyvel,vThres,msMinDur);
    
    % remove the ones that occurr closely together (overshoot)
    right_numSacs = size(right_sacRaw,1);
    minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
    right_interSac = right_sacRaw(2:end,1)- right_sacRaw(1:end-1,2);
    right_sac = right_sacRaw([1; find(right_interSac > minInterSamples)+1],:);
    fprintf('%d rejected for close spacing\n', right_numSacs - size(right_sac,1));
    
    fprintf('%d saccades detected\n', size(right_sac,1));
    
    % saved detected saccades into s
    s.sacsRaw = right_sacRaw;
    s.sacs = right_sac;
    s.sacDetectRadius  = right_radius;
    s.eyeInfo.vThres   = vThres;
    s.eyeInfo.msMinDur = msMinDur;
    s.time = eye_ts_on_right_epoched;
    s.xyPos = right_xypos;
    
    fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
        size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));
    
    % look at detected saccades
    figure(6); clf; set(gcf,'Color', 'w');
    msSacStats1(s);
    title('Angular distribution RIGHT');
    
    
    
    
    
    
    
    
    
    %% Get microsaccades
    % Both
    [both_sacRaw,both_radius] = microsacc(both_xypos,both_xyvel,vThres,msMinDur);
    % Right
    [right_sacRaw,right_radius] = microsacc(right_xypos,right_xyvel,vThres,msMinDur);
    % Left
    [left_sacRaw,left_radius] = microsacc(left_xypos,left_xyvel,vThres,msMinDur);
    % Off
    [off_sacRaw,off_radius] = microsacc(off_xypos,off_xyvel,vThres,msMinDur);
    
    
    
    %% Remove the ones that occurr closely together (overshoot)
    % Both
    both_numSacs = size(both_sacRaw,1);
    both_interSac = both_sacRaw(2:end,1)- both_sacRaw(1:end-1,2);
    both_sac = both_sacRaw([1; find(both_interSac > minInterSamples)+1],:);
    
    % Right
    right_numSacs = size(right_sacRaw,1);
    right_interSac = right_sacRaw(2:end,1)- right_sacRaw(1:end-1,2);
    right_sac = right_sacRaw([1; find(right_interSac > minInterSamples)+1],:);
    
    % Left
    left_numSacs = size(left_sacRaw,1);
    left_interSac = left_sacRaw(2:end,1)- left_sacRaw(1:end-1,2);
    left_sac = left_sacRaw([1; find(left_interSac > minInterSamples)+1],:);
    
    % Left
    off_numSacs = size(off_sacRaw,1);
    off_interSac = off_sacRaw(2:end,1)- off_sacRaw(1:end-1,2);
    off_sac = off_sacRaw([1; find(off_interSac > minInterSamples)+1],:);
    
    %% Prepare to plot sample and saccade distribution (95% confidence)
    
    % ALL trials
    all_traceMean = nanmean(s.xyPos);
    all_traceMedian = nanmedian(s.xyPos);
    all_traceCov = cov(s.xyPos(~isnan(s.xyPos(:,1)),:));
    
    % Both
    both_traceMean = nanmean(both_xypos);
    both_traceMedian = nanmedian(both_xypos);
    both_traceCov = cov(both_xypos(~isnan(both_xypos(:,1)),:));
    
    % RIGHT trials
    right_traceMean = nanmean(right_xypos);
    right_traceMedian = nanmedian(right_xypos);
    right_traceCov = cov(right_xypos(~isnan(right_xypos(:,1)),:));
    
    % LEFT trials
    left_traceMean = nanmean(left_xypos);
    left_traceMedian = nanmedian(left_xypos);
    left_traceCov = cov(left_xypos(~isnan(left_xypos(:,1)),:));
    
    % OFF trials
    off_traceMean = nanmean(off_xypos);
    off_traceMedian = nanmedian(off_xypos);
    off_traceCov = cov(off_xypos(~isnan(off_xypos(:,1)),:));
    
    
    
    %% Plot all samples
    
    figure(7); clf; set(gcf,'Color', 'w');
    subplot(1,2,1); hold on;
    
    % distribution of samples
    plot(s.xyPos(:,1),s.xyPos(:,2),'.','markersize',1);
    % 95% confidence ellipse
    error_ellipse(all_traceCov,all_traceMean,'conf',0.95,'color','k');
    % median
    plot(all_traceMedian(1),all_traceMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Eye position of all samples');
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    
    subplot(1,2,2); hold on; % distribution of saccades
    dxy = s.sacs(:,6:7);
    sacMean = mean(dxy);
    sacMedian = median(dxy);
    sacCov = cov(dxy);
    % plot all samples
    plot(dxy(:,1),dxy(:,2),'.','markersize',2);
    % 95% confidence ellipse
    error_ellipse(sacCov,sacMean,'conf',0.95,'color','k');
    % median
    plot(sacMedian(1),sacMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Saccade vectors');
    
    
    %% Plot off baseline samples
    
    figure(8); clf; set(gcf,'Color', 'w');
    subplot(1,2,1); hold on;
    
    % distribution of samples
    plot(off_xypos(:,1),off_xypos(:,2),'.','markersize',1);
    % 95% confidence ellipse
    error_ellipse(off_traceCov,off_traceMean,'conf',0.95,'color','k');
    % median
    plot(off_traceMedian(1),off_traceMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Eye position of OFF samples');
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    
    subplot(1,2,2); hold on; % distribution of saccades
    dxy = off_sac(:,6:7);
    sacMean = mean(dxy);
    sacMedian = median(dxy);
    sacCov = cov(dxy);
    % plot all samples
    plot(dxy(:,1),dxy(:,2),'.','markersize',2);
    % 95% confidence ellipse
    error_ellipse(sacCov,sacMean,'conf',0.95,'color','k');
    % median
    plot(sacMedian(1),sacMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Saccade vectors');
    
    
    
    
    %% Plot Both samples
    figure(9); clf; set(gcf,'Color', 'w');
    subplot(1,2,1); hold on; % distribution of samples
    
    plot(both_xypos(:,1),both_xypos(:,2),'.','markersize',1);
    % 95% confidence ellipse
    error_ellipse(both_traceCov,both_traceMean,'conf',0.95,'color','k');
    % median
    plot(both_traceMedian(1),both_traceMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Eye position of Attention Both samples ');
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    
    subplot(1,2,2); hold on; % distribution of saccades
    dxy = both_sac(:,6:7);
    sacMean = mean(dxy);
    sacMedian = median(dxy);
    sacCov = cov(dxy);
    % plot all samples
    plot(dxy(:,1),dxy(:,2),'.','markersize',2);
    % 95% confidence ellipse
    error_ellipse(sacCov,sacMean,'conf',0.95,'color','k');
    % median
    plot(sacMedian(1),sacMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Saccade vectors');
    
    
    %% Plot Right samples
    figure(10); clf; set(gcf,'Color', 'w');
    subplot(1,2,1); hold on; % distribution of samples
    
    plot(right_xypos(:,1),right_xypos(:,2),'.','markersize',1);
    % 95% confidence ellipse
    error_ellipse(right_traceCov,right_traceMean,'conf',0.95,'color','k');
    % median
    plot(right_traceMedian(1),right_traceMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Eye position of Attention Right samples');
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    
    subplot(1,2,2); hold on; % distribution of saccades
    dxy = right_sac(:,6:7);
    sacMean = mean(dxy);
    sacMedian = median(dxy);
    sacCov = cov(dxy);
    % plot all samples
    plot(dxy(:,1),dxy(:,2),'.','markersize',2);
    % 95% confidence ellipse
    error_ellipse(sacCov,sacMean,'conf',0.95,'color','k');
    % median
    plot(sacMedian(1),sacMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Saccade vectors');
    
    
    %% Plot Left samples
    figure(11); clf; set(gcf,'Color', 'w');
    subplot(1,2,1); hold on; % distribution of samples
    
    plot(left_xypos(:,1),left_xypos(:,2),'.','markersize',1);
    % 95% confidence ellipse
    error_ellipse(left_traceCov,left_traceMean,'conf',0.95,'color','k');
    % median
    plot(left_traceMedian(1),left_traceMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Eye position of Attention Left samples');
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    
    subplot(1,2,2); hold on; % distribution of saccades
    dxy = left_sac(:,6:7);
    sacMean = mean(dxy);
    sacMedian = median(dxy);
    sacCov = cov(dxy);
    % plot all samples
    plot(dxy(:,1),dxy(:,2),'.','markersize',2);
    % 95% confidence ellipse
    error_ellipse(sacCov,sacMean,'conf',0.95,'color','k');
    % median
    plot(sacMedian(1),sacMedian(2),'+r','markersize',5);
    grid on; axis square;
    xlabel('Horizontal (deg)'); ylabel('Vertical (deg)');
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Saccade vectors');
    
    
    
    %% Plot distributions of saccades on top of eachother
    
    % Both
    both_dxy = both_sac(:,6:7);
    both_sacMean = mean(both_dxy);
    both_sacMedian = median(both_dxy);
    both_sacCov = cov(both_dxy);
    
    % Left
    left_dxy = left_sac(:,6:7);
    left_sacMean = mean(left_dxy);
    left_sacMedian = median(left_dxy);
    left_sacCov = cov(left_dxy);
    
    % Right
    right_dxy = right_sac(:,6:7);
    right_sacMean = mean(right_dxy);
    right_sacMedian = median(right_dxy);
    right_sacCov = cov(right_dxy);
    
    % Off
    off_dxy = off_sac(:,6:7);
    off_sacMean = mean(off_dxy);
    off_sacMedian = median(off_dxy);
    off_sacCov = cov(off_dxy);
    
    
    
    figure(12); clf; set(gcf,'Color', 'w');
    
    
    subplot(1,2,1); hold on; % distribution of samples
    
    % 95% confidence ellipse
    error_ellipse(both_traceCov,both_traceMean,'conf',0.95,'color','r'); hold on;
    % error_ellipse(left_traceCov,left_traceMean,'conf',0.95,'color','g');
    % error_ellipse(right_traceCov,right_traceMean,'conf',0.95,'color','b');
    error_ellipse(off_traceCov,off_traceMean,'conf',0.95,'color','k');
    %
    % median
    plot(both_traceMedian(1),both_traceMedian(2),'+r','markersize',5); hold on;
    % plot(left_traceMedian(1),left_traceMedian(2),'+g','markersize',5);
    % plot(right_traceMedian(1),right_traceMedian(2),'+b','markersize',5);
    plot(off_traceMedian(1),off_traceMedian(2),'+k','markersize',5);
    
    
    grid on; axis square;
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Eye position', 'FontSize',18);
    set(gca,  'FontSize',18); xlabel('Horizontal (deg)'); set(gca,  'FontSize',18); ylabel('Vertical (deg)');
    legend('Both', 'Blank'); % ,
    
    
    
    subplot(1,2,2); hold on;
    
    % 95% confidence ellipse
    error_ellipse(both_sacCov,both_sacMean,'conf',0.95,'color','r'); hold on;
    % error_ellipse(left_sacCov,left_sacMean,'conf',0.95,'color','g');
    % error_ellipse(right_sacCov,right_sacMean,'conf',0.95,'color','b');
    error_ellipse(off_sacCov,off_sacMean,'conf',0.95,'color','k');
    
    % plot both samples
    plot(both_dxy(:,1),both_dxy(:,2),'r.','markersize',2); hold on;
    % plot(left_dxy(:,1),left_dxy(:,2),'g.','markersize',2);
    % plot(right_dxy(:,1),right_dxy(:,2),'b.','markersize',2);
    plot(off_dxy(:,1),off_dxy(:,2), 'k.','markersize',2);
    
    
    % median
    plot(both_sacMedian(1),both_sacMedian(2),'+r','markersize',5); hold on;
    % plot(left_sacMedian(1),left_sacMedian(2),'+g','markersize',5);
    % plot(right_sacMedian(1),right_sacMedian(2),'+b','markersize',5);
    plot(off_sacMedian(1),off_sacMedian(2),'+k','markersize',5);
    
    
    grid on; axis square;
    set(gca,  'FontSize',18); xlabel('Horizontal (deg)'); set(gca,  'FontSize',18); ylabel('Vertical (deg)');
    xlim(5*[-1,1]); ylim(5*[-1,1]);
    title('Saccade vectors' ,  'FontSize',18);
    legend('Both', 'Blank'); %'Left', 'Right',
    
    
    
    
    
    
    
    
    
    %% Export images
    if save_images
        if delete_first_last
            cd(datapath);
            hgexport(1,fullfile(datapath,'eyetracking_fig1_xypositions_no_firstlast.eps'));
            hgexport(2,fullfile(datapath,'eyetracking_fig2_microsaccades_all_no_firstlast.eps'));
            hgexport(3,fullfile(datapath,'eyetracking_fig3_microsaccades_both_no_firstlast.eps'));
            hgexport(4,fullfile(datapath,'eyetracking_fig4_microsaccades_off_no_firstlast.eps'));
            hgexport(5,fullfile(datapath,'eyetracking_fig5_microsaccades_left_no_firstlast.eps'));
            hgexport(6,fullfile(datapath,'eyetracking_fig6_microsaccades_right_no_firstlast.eps'));
            hgexport(7,fullfile(datapath,'eyetracking_fig7_xypos_ms_all_no_firstlast.eps'));
            hgexport(8,fullfile(datapath,'eyetracking_fig8_xypos_ms_off_no_firstlast.eps'));
            hgexport(9,fullfile(datapath,'eyetracking_fig9_xypos_ms_both_no_firstlast.eps'));
            hgexport(10,fullfile(datapath,'eyetracking_fig10_xypos_ms_right_no_firstlast.eps'));
            hgexport(11,fullfile(datapath,'eyetracking_fig11_xypos_ms_left_no_firstlast.eps'));
            hgexport(12,fullfile(datapath,'eyetracking_fig12_xypos_ms_combined_no_firstlast.eps'));
            
        else
            
            hgexport(1,fullfile(datapath,'eyetracking_fig1_xypositions.eps'));
            hgexport(2,fullfile(datapath,'eyetracking_fig2_microsaccades_all.eps'));
            hgexport(3,fullfile(datapath,'eyetracking_fig3_microsaccades_both.eps'));
            hgexport(4,fullfile(datapath,'eyetracking_fig4_microsaccades_off.eps'));
            hgexport(5,fullfile(datapath,'eyetracking_fig5_microsaccades_left.eps'));
            hgexport(6,fullfile(datapath,'eyetracking_fig6_microsaccades_right.eps'));
            hgexport(7,fullfile(datapath,'eyetracking_fig7_xypos_ms_all.eps'));
            hgexport(8,fullfile(datapath,'eyetracking_fig8_xypos_ms_off.eps'));
            hgexport(9,fullfile(datapath,'eyetracking_fig9_xypos_ms_both.eps'));
            hgexport(10,fullfile(datapath,'eyetracking_fig10_xypos_ms_right.eps'));
            hgexport(11,fullfile(datapath,'eyetracking_fig11_xypos_ms_left.eps'));
            hgexport(12,fullfile(datapath,'eyetracking_fig12_xypos_ms_combined.eps'));
        end
    end
    
end
