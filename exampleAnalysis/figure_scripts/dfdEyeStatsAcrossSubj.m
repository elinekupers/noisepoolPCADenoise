%% dfdEyeStatsAcrossSubj

% Script to load saccade statistics saved by dfdEyeScript and make bargraph
% of stats across subjects

subjects             = 6:8;
dataDir              = fullfile(dfdRootPath, 'exampleAnalysis', 'data');   % Where to save data?
figureDir            = fullfile(dfdRootPath, 'exampleAnalysis', 'figures_rm1epoch');% Where to save images?
saveFigures          = true;   % Save figures in the figure folder?
dataAll              = [];

%% Loops over datasets
% cd(datapath)

for whichSubject = 1:length(subjects);
    dataAll{whichSubject} = load(fullfile(dataDir,sprintf('s0%d_freq_mdAmpl.mat',subjects(whichSubject))));   
end

% Take mean freq and sd
bothFreq = [];
LeftFreq = [];
RightFreq = [];
BlankFreq = [];

bothAmplX = [];
bothAmplY = [];

leftAmplX = [];
leftAmplY = [];

rightAmplX = [];
rightAmplY = [];

blankAmplX = [];
blankAmplY = [];

bothMS = [];
leftMS    = [];
rightMS   = [];


dataAll =dataAll(~cellfun('isempty',dataAll)); 

for subject_num = [1:3];
    
    bothFreq(subject_num,:) = dataAll{subject_num}.stats.freqBoth;
    LeftFreq(subject_num,:) = dataAll{subject_num}.stats.freqLeft;
    RightFreq(subject_num,:) = dataAll{subject_num}.stats.freqRight;

    blankFreq(subject_num,:) = dataAll{subject_num}.stats.freqBlank;
    
    dataAll{1}.stats.allData

    
    bothAmplX(subject_num,:) = dataAll{subject_num}.stats(1).mdAmplBoth(1);
    bothAmplY(subject_num,:) = dataAll{subject_num}.stats(1).mdAmplBoth(2);
    leftAmplX(subject_num,:) = dataAll{subject_num}.stats(2).mdAmplLeft(1);
    leftAmplY(subject_num,:) = dataAll{subject_num}.stats(2).mdAmplLeft(2);
    rightAmplX(subject_num,:) = dataAll{subject_num}.stats(3).mdAmplRight(1);
    rightAmplY(subject_num,:) = dataAll{subject_num}.stats(3).mdAmplRight(2);
    blankAmplX(subject_num,:) = dataAll{subject_num}.stats(4).mdAmplBlank(1);
    blankAmplY(subject_num,:) = dataAll{subject_num}.stats(4).mdAmplBlank(2);
    

end

% Take mean xpos, ypos and sd xpos, ypos
mnB = mean(bothFreq);
sdB = std(bothFreq)/sqrt(length(3));

mnL = mean(leftFreq);
sdL = std(leftFreq)/sqrt(length(3));

mnR = mean(rightFreq);
sdR = std(rightFreq)/sqrt(length(3));

mnBl = mean(blanksFreq);
sdBl = std(blanksFreq)/sqrt(length(3));

mnBX = mean(bothAmplX);
sdBX = std(bothAmplX)/sqrt(length(3));
mnBY = mean(bothAmplY);
sdBY = std(bothAmplY)/sqrt(length(3));

mnLX = mean(leftAmplX);
sdLX = std(leftAmplX)/sqrt(length(3));
mnLY = mean(leftAmplY);
sdLY = std(leftAmplY)/sqrt(length(3));

mnRX = mean(rightAmplX);
sdRX = std(rightAmplX)/sqrt(length(3));
mnRY = mean(rightAmplY);
sdRY = std(rightAmplY)/sqrt(length(3));

mnBLX = mean(blankAmplX);
sdBLX = std(blankAmplX)/sqrt(length(3));
mnBLY = mean(blankAmplY);
sdBLY = std(blankAmplY)/sqrt(length(3));

figure; clf; 
subplot(1,1,1);
bar_width = .5;
% mean and sem across subjects
bar([1,2,3], [mnB,mnL,mnR],bar_width , 'EdgeColor','none'); hold on;
errorbar2([1,2,3], [mnB,mnL,mnR], [sdB,sdL,sdR],1,'k-');
set(gca,'XTickLabel',{'Both', 'Left', 'Right'})
% format figure and make things pretty
set(gca,'xlim',[0.2,3+0.8],'ylim',[0,700]);
makeprettyaxes(gca,9,9);
ylabel('Frequency','FontSize',18)
xlabel('Conditions','FontSize',18)

% hgexport(gcf,fullfile(datapath,'across_subjects_eyetracking_stats_bar.eps'));

rgb_purple    = [93,12,139]/255;
satValues       = 1-linspace(0.1,1,9);
colorRGB_pink   = varysat(rgb_purple,satValues);
colorRGB_pink   = squeeze(colorRGB_pink);
colorRGB_pink   = colorRGB_pink';

[~, idx] = sort(bothFreq);
bothFreq = bothFreq(idx);
LeftFreq   = LeftFreq(idx);
blanks_freq  = blanks_freq(idx);


% x = [ones(1,size(grating_freq,1)),2*ones(1,size(noise_freq,1)),3*ones(1,size(blanks_freq,1))];
y = [bothFreq,LeftFreq,blanks_freq];

cmap = colorRGB_pink;


figure; clf; 
subplot(1,1,1);
for k = 1:size(y,1)
    plot([1:1:3], [y(k,1),y(k,2),y(k,3)], 'o-','color', cmap(k,:), 'linewidth',4); hold on;
end

set(gca,'XTickLabel',{'','Gratings', '','Noise','', 'Blank',''}, 'FontSize',20)
% format figure and make things pretty
set(gca,'xlim',[0.2,3+0.8],'ylim',[0,1200]);
makeprettyaxes(gca,9,9);
ylabel('Frequency','FontSize',20)
xlabel('Conditions','FontSize',20)

hgexport(gcf,fullfile(datapath,'across_subjects_eyetracking_stats_scatter.eps'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Sanity check plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

theseSacs = s.sacs;
[rho,theta] = msGetThetaRho(s.xyPos,theseSacs);
peakVel = theseSacs(:,3);

% plot main sequence
subplot(2,2,1);
loglog(rho,peakVel,'or','markersize',2);
axis square, grid on
xlabel('Amplitude')
ylabel('Peak velocity')
title('Main Sequence')

% plot intersaccade interval
subplot(2,2,2);
sacTimes = s.time(theseSacs(:,1))/1000;
sacTimesInter = diff(sacTimes);

hist(sacTimesInter,50);
axis square
xlabel('IMSI duration (sec)')
ylabel('Frequency')
title('Inter-microsaccade interval distribution')

% plot individual microsaccade vectors
subplot(2,2,3);
h=polar2(theta,rho,'o');
set(h,'markersize',2);
title('Microsaccade vectors')

% circular distribution
subplot(2,2,4);
rose(theta,(0:10:360)/180*pi);
title('Angular distribution')




