function dfdEyeScript(subjects)

% 
% dfdEyeScript(subjects)
%
% INPUTS:
% subjects    : Number of subjects one would like to denoise
%
% DESCRIPTION: Function to analyze eyetracking data recorded during MEG visual steady
% state experiment. 
%
% AUTHORS. YEAR. TITLE. JOURNAL.

% ------------------------------------------------------------------------
% ---------------- Check options and load in toolboxes -------------------
% ------------------------------------------------------------------------

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

% ------------------------------------------------------------------------
% ------------------ Define paths and load in data -----------------------
% ------------------------------------------------------------------------
for whichSubject = subjects
edffile = dir(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye','s0%d_eyelink.edf'),whichSubject));

if saveEyd
    eyd = mglEyelinkEDFRead(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye',edffile.name)));
    save(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye', 's0%d_eyd.mat'),whichSubject));
else
    eyd = load(sprintf(fullfile(dfdRootPath, 'exampleAnalysis', 'data', 'eye','s0*_eyd.mat'),whichSubject));
end

% ------------------------------------------------------------------------
% --------------- Get eye traces (blinks are removed) --------------------
% ------------------------------------------------------------------------

% Figure out the start and end time of your experiment, which will be a
% an argument for the msGetEyeData function. 
% Blink window is pretty conservative right now. could be [-0.1,0.1] (100
% ms either way)
startTime   = 8;
endTime     = 0;
timeLims    = [eyd.messages(startTime).time(1), eyd.messages(end-endTime).time];
blinkWinSec = [0.2 0.35];
s = msGetEyeData(eyd,timeLims,blinkWinSec);

% ------------------------------------------------------------------------
% ------------------- Get triggers and timings ---------------------------
% ------------------------------------------------------------------------

% Use MEG messages to get the eyemovements for separate conditions

% ---------------------- Delete irrelevant messages ----------------------
eyd.messages = eyd.messages(1,startTime:end);

% ------------- Make a matrix for trigger nr and time stamp --------------
triggers = zeros(size(eyd.messages(1,:),2),2);

% ------------------ Get trigger nr and time stamp -----------------------
for ii = 1:size(eyd.messages(1,:),2); % Last triggers while quiting the program are getting deleted
    triggers(ii,1) = str2num(eyd.messages(1,ii).message(14)); % Trigger nr
    triggers(ii,2) = eyd.messages(1,ii).time(1)-s.timeRaw(1); % Get timing and set time to zero
end

% ------------------ Add triggers for the blank periods ------------------

onsets = ssmeg_trigger_2_onsets(triggers, whichSubject, 'eye');

% Try to use "onsets = ssmeg_trigger_2_onsets(triggers, which_subject)" to
% get also triggers for blank periods.
% then use meg_make_epochs with s.time as raw ts 


numTimePoints = 1000;
eye_ts = s.time(triggers(1:12:end,2) + (0:numTimePoints-1));

% ------------------ Make design matrix ----------------------------------
designMatrix = zeros(size(eyd.messages(1,:),2),3);
designMatrix(conditions==1,1) = 1; % condition 1 is full field
designMatrix(conditions==5,2) = 1; % condition 5 is left field
designMatrix(conditions==7,3) = 1; % condition 7 is right field







% 
% 


return

%% Make 1s epoch (Chunk triggers)
cycles_per_s = 6; % six triggers in one second, so 36 in one 6 second trial

% Get start value of every 1-s 'on' epoch
epoch_start_both_att_on_ind  = trig_both_att_ind(1:cycles_per_s:end); 
epoch_start_right_att_on_ind = trig_right_att_ind(1:cycles_per_s:end);
epoch_start_left_att_on_ind  = trig_left_att_ind(1:cycles_per_s:end);

% Check number of epochs
assert(size(epoch_start_both_att_on_ind,1) + size(epoch_start_right_att_on_ind,1)+size(epoch_start_left_att_on_ind,1)==6*size(trig_probe_onset_ind,1))

num_epoch_time_pts = mode(diff(epoch_start_both_att_on_ind)); % frame duration should be the same for every condition

% Get start value of every 1-s 'off' epoch
epoch_start_both_att_off_ind  = epoch_start_both_att_on_ind - num_epoch_time_pts;
epoch_start_right_att_off_ind = epoch_start_right_att_on_ind -  num_epoch_time_pts;
epoch_start_left_att_off_ind  = epoch_start_left_att_on_ind - num_epoch_time_pts;

%% Define in epochs in the timeseries

ep_both = length(epoch_start_both_att_on_ind);
ep_right = length(epoch_start_right_att_on_ind);
ep_left = length(epoch_start_left_att_on_ind);

% Make new arrays
eye_ts_on_both_epoched   = zeros(num_epoch_time_pts, ep_both); % Time x Epochs
eye_ts_on_right_epoched  = zeros(num_epoch_time_pts, ep_right); 
eye_ts_on_left_epoched   = zeros(num_epoch_time_pts, ep_left);
eye_ts_off_both_epoched  = zeros(num_epoch_time_pts, ep_both); 
eye_ts_off_right_epoched = zeros(num_epoch_time_pts, ep_right); 
eye_ts_off_left_epoched  = zeros(num_epoch_time_pts, ep_left);


%% Loop over channels and epochs of every conditiom

% Both
for epoch = 1:ep_both
    
    % Time 'on'
    eye_ts_on_both_epoched(:, epoch) = s.time(epoch_start_both_att_on_ind(epoch) + (0:num_epoch_time_pts-1));    
    % Time 'off'
    eye_ts_off_both_epoched(:, epoch) = s.time(epoch_start_both_att_off_ind(epoch) + (0:num_epoch_time_pts-1));
    
    % Position 'on'
    eye_xpos_on_both_epoched(:, epoch) = s.xyPos(epoch_start_both_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_ypos_on_both_epoched(:, epoch) = s.xyPos(epoch_start_both_att_on_ind(epoch) + (0:num_epoch_time_pts-1),2);
    % Position 'off' 
    eye_xpos_off_both_epoched(:, epoch) = s.xyPos(epoch_start_both_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_ypos_off_both_epoched(:, epoch) = s.xyPos(epoch_start_both_att_off_ind(epoch) + (0:num_epoch_time_pts-1),2);
    
    % Velocity on
    eye_xvel_on_both_epoched(:, epoch) = s.xyVel(epoch_start_both_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_yvel_on_both_epoched(:, epoch) = s.xyVel(epoch_start_both_att_on_ind(epoch) + (0:num_epoch_time_pts-1),2);
    % Velocity off
    eye_xvel_off_both_epoched(:, epoch) = s.xyVel(epoch_start_both_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_yvel_off_both_epoched(:, epoch) = s.xyVel(epoch_start_both_att_off_ind(epoch) + (0:num_epoch_time_pts-1),2);
end

% RIGHT
for epoch = 1:ep_right
    
    % Time on
    eye_ts_on_right_epoched(:, epoch) = s.time(epoch_start_right_att_on_ind(epoch) + (0:num_epoch_time_pts-1));
    % Time off
    eye_ts_off_right_epoched(:, epoch) = s.time(epoch_start_right_att_off_ind(epoch) + (0:num_epoch_time_pts-1));
    
    % Position on
    eye_xpos_on_right_epoched(:, epoch) = s.xyPos(epoch_start_right_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_ypos_on_right_epoched(:, epoch) = s.xyPos(epoch_start_right_att_on_ind(epoch) + (0:num_epoch_time_pts-1),2);
    % Position off
    eye_xpos_off_right_epoched(:, epoch) = s.xyPos(epoch_start_right_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_ypos_off_right_epoched(:, epoch) = s.xyPos(epoch_start_right_att_off_ind(epoch) + (0:num_epoch_time_pts-1),2);
    
    % Velocity on
    eye_xvel_on_right_epoched(:, epoch) = s.xyVel(epoch_start_right_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_yvel_on_right_epoched(:, epoch) = s.xyVel(epoch_start_right_att_on_ind(epoch) + (0:num_epoch_time_pts-1),2);
    % Velocity off
    eye_xvel_off_right_epoched(:, epoch) = s.xyVel(epoch_start_right_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_yvel_off_right_epoched(:, epoch) = s.xyVel(epoch_start_right_att_off_ind(epoch) + (0:num_epoch_time_pts-1),2); 
end

% LEFT
for epoch = 1:ep_left
    
    % Time on
    eye_ts_on_left_epoched(:, epoch) = s.time(epoch_start_left_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    % Time off
    eye_ts_off_left_epoched(:, epoch) = s.time(epoch_start_left_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    
    % Position on
    eye_xpos_on_left_epoched(:, epoch) = s.xyPos(epoch_start_left_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_ypos_on_left_epoched(:, epoch) = s.xyPos(epoch_start_left_att_on_ind(epoch) + (0:num_epoch_time_pts-1),2);
    % Position off
    eye_xpos_off_left_epoched(:, epoch) = s.xyPos(epoch_start_left_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_ypos_off_left_epoched(:, epoch) = s.xyPos(epoch_start_left_att_off_ind(epoch) + (0:num_epoch_time_pts-1),2);
    
    % Velocity on
    eye_xvel_on_left_epoched(:, epoch) = s.xyVel(epoch_start_left_att_on_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_yvel_on_left_epoched(:, epoch) = s.xyVel(epoch_start_left_att_on_ind(epoch) + (0:num_epoch_time_pts-1),2);
    % Velocity off
    eye_xvel_off_left_epoched(:, epoch) = s.xyVel(epoch_start_left_att_off_ind(epoch) + (0:num_epoch_time_pts-1),1);
    eye_yvel_off_left_epoched(:, epoch) = s.xyVel(epoch_start_left_att_off_ind(epoch) + (0:num_epoch_time_pts-1),2);

end

% Concatenate the off time series
eye_ts_off_epoched = cat(2, eye_ts_off_both_epoched,eye_ts_off_right_epoched,eye_ts_off_left_epoched);
eye_xpos_off_epoched = cat(2, eye_xpos_off_both_epoched,eye_xpos_off_right_epoched,eye_xpos_off_left_epoched);
eye_ypos_off_epoched = cat(2, eye_ypos_off_both_epoched,eye_ypos_off_right_epoched,eye_ypos_off_left_epoched);
eye_xvel_off_epoched = cat(2, eye_xvel_off_both_epoched,eye_xvel_off_right_epoched,eye_xvel_off_left_epoched);
eye_yvel_off_epoched = cat(2, eye_yvel_off_both_epoched,eye_yvel_off_right_epoched,eye_yvel_off_left_epoched);


%% Delete first and last epoch?

if delete_first_last
    
    % We have a different amount of epochs per condition, so we now simply
    % devide the keep_epochs per condition.
    
    keep_epochs_both = true(1,ep_both);
    keep_epochs_both(1:6:end) = false;
    keep_epochs_both(6:6:end) = false;
    
    keep_epochs_left = true(1,ep_left);
    keep_epochs_left(1:6:end) = false;
    keep_epochs_left(6:6:end) = false;
    
    keep_epochs_right = true(1,ep_right);
    keep_epochs_right(1:6:end) = false;
    keep_epochs_right(6:6:end) = false;
    
    ep_off = size(eye_ts_off_epoched,2);
    keep_epochs_off = true(1,ep_off);
    keep_epochs_off(1:6:end) = false;
    keep_epochs_off(6:6:end) = false;
    
    
    % Both
    eye_ts_on_both_epoched   = eye_ts_on_both_epoched(:, keep_epochs_both);    
    eye_xpos_on_both_epoched = eye_xpos_on_both_epoched(:, keep_epochs_both);
    eye_ypos_on_both_epoched = eye_ypos_on_both_epoched(:, keep_epochs_both);
    eye_xvel_on_both_epoched = eye_xvel_on_both_epoched(:, keep_epochs_both);
    eye_yvel_on_both_epoched = eye_yvel_on_both_epoched(:, keep_epochs_both);
   
    % Left
    eye_ts_on_left_epoched   = eye_ts_on_left_epoched(:,keep_epochs_left);    
    eye_xpos_on_left_epoched = eye_xpos_on_left_epoched(:,keep_epochs_left);
    eye_ypos_on_left_epoched = eye_ypos_on_left_epoched(:,keep_epochs_left);
    eye_xvel_on_left_epoched = eye_xvel_on_left_epoched(:,keep_epochs_left);
    eye_yvel_on_left_epoched = eye_yvel_on_left_epoched(:,keep_epochs_left);
    
    % Right
    eye_ts_on_right_epoched   = eye_ts_on_right_epoched(:,keep_epochs_right);    
    eye_xpos_on_right_epoched = eye_xpos_on_right_epoched(:,keep_epochs_right);
    eye_ypos_on_right_epoched = eye_ypos_on_right_epoched(:,keep_epochs_right);
    eye_xvel_on_right_epoched = eye_xvel_on_right_epoched(:,keep_epochs_right);
    eye_yvel_on_right_epoched = eye_yvel_on_right_epoched(:,keep_epochs_right);
    
    % Off
    eye_ts_off_epoched   = eye_ts_off_epoched(:,keep_epochs_off);
    eye_xpos_off_epoched = eye_xpos_off_epoched(:,keep_epochs_off);
    eye_ypos_off_epoched = eye_ypos_off_epoched(:,keep_epochs_off);
    eye_xvel_off_epoched = eye_xvel_off_epoched(:,keep_epochs_off);
    eye_yvel_off_epoched = eye_yvel_off_epoched(:,keep_epochs_off);  
end


%% Save data?
if save_data
    
    cd(datapath)
    data_array = {eye_ts_on_both_epoched,eye_xpos_on_both_epoched,eye_ypos_on_both_epoched,eye_xvel_on_both_epoched,eye_yvel_on_both_epoched,...
        eye_ts_on_left_epoched, eye_xpos_on_left_epoched, eye_ypos_on_left_epoched, eye_xvel_on_left_epoched, eye_yvel_on_left_epoched, ...
        eye_ts_on_right_epoched, eye_xpos_on_right_epoched, eye_ypos_on_right_epoched, eye_xvel_on_right_epoched, eye_yvel_on_right_epoched, ...
        eye_ts_off_epoched, eye_xpos_off_epoched, eye_ypos_off_epoched, eye_xvel_off_epoched, eye_yvel_off_epoched};
    save(fullfile(datapath,'epoched_eye_data'),'data_array');
     
end




%% ================================================= 
%  ===== Plot eye traces for visual inspection =====
%  ================================================= 

deglims = 20;
figure(1); clf; set(gcf,'Color', 'w');

%% X COORDINATES
subplot(2,2,1);

% All eyetracking data
plot(s.time/1000,s.xyPos(:,1), 'k'); hold on;

% All 'off' baselines
plot(eye_ts_off_epoched/1000,eye_xpos_off_epoched, 'Color', [.8 .8 .8]);

% All 'on' flickers
plot(eye_ts_on_both_epoched/1000,eye_xpos_on_both_epoched, 'r'); 
plot(eye_ts_on_left_epoched/1000,eye_xpos_on_left_epoched, 'b'); 
plot(eye_ts_on_right_epoched/1000,eye_xpos_on_right_epoched, 'g'); 


grid on;
ylim(deglims*[-1,1]); xlabel('Time (s)'); ylabel('X (deg)');


%% Y COORDINATES

subplot(2,2,3);

% All eyetracking data
plot(s.time/1000,s.xyPos(:,2),'k'); hold on;

% All 'off' baselines
plot(eye_ts_off_epoched/1000,eye_ypos_off_epoched, 'Color', [.8 .8 .8]);

% All 'on' flickers
plot(eye_ts_on_both_epoched/1000,eye_ypos_on_both_epoched, 'r');
plot(eye_ts_on_left_epoched/1000,eye_ypos_on_left_epoched, 'b');
plot(eye_ts_on_right_epoched/1000,eye_ypos_on_right_epoched, 'g');


grid on;
ylim(deglims*[-1,1]); xlabel('Time (s)'); ylabel('Y (deg)');

%% XY ON GRID
% subplot(2,2,[2,4]);

figure; clf; set(gcf, 'Color', 'w')
% All eyetracking data
plot(s.xyPos(:,1),s.xyPos(:,2), 'k'); axis square; grid on; hold on;

% All 'off' baselines
% plot(eye_xpos_off_epoched,eye_ypos_off_epoched, 'Color', [.8 .8 .8]);

% All 'on' flickers'
% plot(eye_xpos_on_both_epoched,eye_ypos_on_both_epoched, 'r');
plot(eye_xpos_on_left_epoched,eye_ypos_on_left_epoched, 'b'); 
plot(eye_xpos_on_right_epoched,eye_ypos_on_right_epoched, 'g'); 


xlim(deglims*[-1,1]); ylim(deglims*[-1,1]);
xlabel('X (deg)', 'Fontsize', 20);  ylabel('Y (deg)', 'Fontsize', 20);
set(gca, 'FontSize', 20);






% figure; clf; set(gcf,'Color', 'w'); % LR
% 
% % All eyetracking data
% plot(s.xyPos(:,1),s.xyPos(:,2), 'Color', [.5 .5 .5]); axis square; grid on; hold on;
% 
% % All 'on' flickers'
% plot(eye_xpos_on_left_epoched,eye_ypos_on_left_epoched, 'b'); 
% plot(eye_xpos_on_right_epoched,eye_ypos_on_right_epoched, 'g'); 
% 
% 
% xlim(deglims*[-1,1]); ylim(deglims*[-1,1]);
% set(gca, 'FontSize', 18); xlabel('X (deg)', 'Fontsize', 18);  set(gca, 'FontSize', 18); ylabel('Y (deg)', 'Fontsize', 18);
% 
% 
% figure; clf; set(gcf,'Color', 'w'); % B-Off
% 
% % All eyetracking data
% plot(s.xyPos(:,1),s.xyPos(:,2), 'Color', [.5 .5 .5]); axis square; grid on; hold on;
% 
% % All 'off' baselines
% plot(eye_xpos_off_epoched,eye_ypos_off_epoched, 'k');
% 
% % All 'on' flickers'
% plot(eye_xpos_on_both_epoched,eye_ypos_on_both_epoched, 'r');
% 
% 
% xlim(deglims*[-1,1]); ylim(deglims*[-1,1]);
% set(gca, 'FontSize', 18); xlabel('X (deg)', 'Fontsize', 18);  set(gca, 'FontSize', 18); ylabel('Y (deg)', 'Fontsize', 18);
% 





%% Concatenate all conditions for XY position and XY velocity
% BOTH
new_array_sz = size(eye_xvel_on_both_epoched,1)*size(eye_xvel_on_both_epoched,2);
both_xyvel = [reshape(eye_xvel_on_both_epoched,1,new_array_sz);reshape(eye_yvel_on_both_epoched,1,new_array_sz)]';

new_array_sz = size(eye_xpos_on_both_epoched,1)*size(eye_xpos_on_both_epoched,2);
both_xypos = [reshape(eye_xpos_on_both_epoched,1,new_array_sz);reshape(eye_ypos_on_both_epoched,1,new_array_sz)]';

% RIGHT 
new_array_sz = size(eye_xvel_on_right_epoched,1)*size(eye_xvel_on_right_epoched,2);
right_xyvel = [reshape(eye_xvel_on_right_epoched,1,new_array_sz);reshape(eye_yvel_on_right_epoched,1,new_array_sz)]';

new_array_sz = size(eye_xpos_on_right_epoched,1)*size(eye_xpos_on_right_epoched,2);
right_xypos = [reshape(eye_xpos_on_right_epoched,1,new_array_sz);reshape(eye_ypos_on_right_epoched,1,new_array_sz)]';

% LEFT
new_array_sz = size(eye_xvel_on_left_epoched,1)*size(eye_xvel_on_left_epoched,2);
left_xyvel = [reshape(eye_xvel_on_left_epoched,1,new_array_sz);reshape(eye_yvel_on_left_epoched,1,new_array_sz)]';

new_array_sz = size(eye_xpos_on_left_epoched,1)*size(eye_xpos_on_left_epoched,2);
left_xypos = [reshape(eye_xpos_on_left_epoched,1,new_array_sz);reshape(eye_ypos_on_left_epoched,1,new_array_sz)]';

% OFF
new_array_sz = size(eye_xvel_off_epoched,1)*size(eye_xvel_off_epoched,2);
off_xyvel = [reshape(eye_xvel_off_epoched,1,new_array_sz);reshape(eye_yvel_off_epoched,1,new_array_sz)]';

new_array_sz = size(eye_xpos_off_epoched,1)*size(eye_xpos_off_epoched,2);
off_xypos = [reshape(eye_xpos_off_epoched,1,new_array_sz);reshape(eye_ypos_off_epoched,1,new_array_sz)]';


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

%% ALL
[sacRaw,radius] = microsacc(s.xyPos,s.xyVel,vThres,msMinDur);

% remove the ones that occurr closely together (overshoot)
numSacs = size(sacRaw,1);
minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
interSac = sacRaw(2:end,1)- sacRaw(1:end-1,2);
sac = sacRaw([1; find(interSac > minInterSamples)+1],:);
fprintf('%d rejected for close spacing\n', numSacs - size(sac,1));

fprintf('%d saccades detected\n', size(sac,1));

% saved detected saccades into s
s.sacsRaw = sacRaw;
s.sacs = sac;
s.sacDetectRadius  = radius;
s.eyeInfo.vThres   = vThres;
s.eyeInfo.msMinDur = msMinDur;

fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
        size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));

% look at detected saccades
figure(2); clf; set(gcf,'Color', 'w');
msSacStats1(s);
title('Angular distribution ALL');

%% BOTH

[both_sacRaw,both_radius] = microsacc(both_xypos,both_xyvel,vThres,msMinDur);

% remove the ones that occurr closely together (overshoot)
both_numSacs = size(both_sacRaw,1);
minInterSamples = ceil(0.01*s.eyeInfo.smpRate);
both_interSac = both_sacRaw(2:end,1)- both_sacRaw(1:end-1,2);
both_sac = both_sacRaw([1; find(both_interSac > minInterSamples)+1],:);
fprintf('%d rejected for close spacing\n', both_numSacs - size(both_sac,1));

fprintf('%d saccades detected\n', size(both_sac,1));

% saved detected saccades into s
s.sacsRaw = both_sacRaw;
s.sacs = both_sac;
s.sacDetectRadius  = both_radius;
s.eyeInfo.vThres   = vThres;
s.eyeInfo.msMinDur = msMinDur;
s.time = eye_ts_on_both_epoched;
s.xyPos = both_xypos;

fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
        size(s.sacs,1), s.sacDetectRadius(1), s.sacDetectRadius(2));

% look at detected saccades
figure(3); clf; set(gcf,'Color', 'w');
msSacStats1(s);
title('Angular distribution BOTH');

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
