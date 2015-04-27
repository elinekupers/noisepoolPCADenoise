function [trl, event] = ssmeg_ft_trial_fun(cfg)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);
event = cfg.event;

% search for "trigger" events
sample = [cfg.event(find(strcmp('trigger', {cfg.event.type}))).sample]';

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
posttrig =  round(cfg.trialdef.post * hdr.Fs);

% look for the combination of a trigger "7" followed by a trigger "64" 
% for each trigger except the last one
trl = [];
for j = 1:(length(value))
    trlbegin = sample(j) + pretrig;
    trlend   = sample(j) + posttrig;
    offset   = pretrig;
    newtrl   = [trlbegin trlend offset];
    trl      = [trl; newtrl];
end