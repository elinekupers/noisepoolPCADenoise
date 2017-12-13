function s = msGetEyeData(edf,timelims,blinkWinSec)
%s = msGetEyeData(edf,varargin)
%gets position and velocity data and saves into s struct
% edf: eyelink matlab struct 
% [timelims]: the start and end times of an experimental block from which
% to pull out eye data (nx2 matrix, where col 1=start times, and col 2=end times)
% parameters:
% blinkWinSec: window in sec for cleaning up blinks. Takes either one value
% (symmetric window) or two values (assymetric window) (Default = [0.2 0.35]

% process input parameters
if ~exist('timelims','var') || isempty(timelims)
     timelims = [edf.messages(1).time(1), edf.messages(end).time];
end
if ~exist('blinkWinSec','var') || isempty(blinkWinSec)
    blinkWinSec = [0.2 0.35];
end

% get sample rate
smpRate = msGetSmpRate(edf);
fprintf('(msGetEyeData) smpRate = %d Hz\n', smpRate);

% filter out blinks from eye data
[xyData,tData] = removeBlinks(edf,blinkWinSec);

for k = 1:size(timelims,1)
    % pull out relevant time points
    tInds  = tData >= timelims(k,1) & tData <= timelims(k,2);
    s(k).timeRaw  = tData(tInds,:);
    s(k).time     = s(k).timeRaw-s(k).timeRaw(1);
    
    s(k).xyPosRaw  = xyData(tInds,:);
    
    % convert to degrees of visual angle
    nSamp  = length(find(tInds));
    s(k).xyPos  = s(k).xyPosRaw - repmat(edf.gazeCoords(3:4)/2, nSamp, 1);
    s(k).xyPos  = s(k).xyPos ./ [edf.gaze.pix2degX(tInds)' edf.gaze.pix2degY(tInds)'];
    % invert y so 0 is bottom left 
    s.xyPos(:,2) = -s.xyPos(:,2);
    %PIX2DEG  = 28.8143575730226;
    %s.xyPos  = s.xyPos ./ PIX2DEG;
    
    % compute velocity
    s(k).xyVel  = vecvel(s(k).xyPos,smpRate,3);
    
    % save eye data information
    s(k).eyeInfo.smpRate     = smpRate;
    s(k).eyeInfo.blinkWinSec = blinkWinSec;
    s(k).eyeInfo.edfFile     = edf.filename;
    s(k).eyeInfo.timeLims    = timelims(k,:);
end
