function [xyData,tData] = removeBlinks(eyd,blinkWinSec)
%[xyData,tData] = removeBlinks(eyd,blinkWinSec)
%remove blinks based on the eye link detection algorithm

% window for taking out blink samples
if ~exist('blinkWinSec','var')
    blinkWinSec = 0.2;
end

% allow for asymmetric window
if length(blinkWinSec) == 1
    blinkWinSec = blinkWinSec*ones(1,2);
end

xyData = [eyd.gaze.x',eyd.gaze.y'];
tData  = eyd.gaze.time';
blinkData = eyd.blinks;

for iBlink = 1:length(blinkData.startTime)
    blinkInd = tData >= blinkData.startTime(iBlink)-blinkWinSec(1)*1000 & ...
               tData <= blinkData.endTime(iBlink)+blinkWinSec(2)*1000;
    xyData(blinkInd,:) = nan;
end