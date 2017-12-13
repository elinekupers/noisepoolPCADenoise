function smpRate = msGetSmpRate(edf)

smpRate = 1000/round(mean(diff(edf.gaze.time(1:10)))); 