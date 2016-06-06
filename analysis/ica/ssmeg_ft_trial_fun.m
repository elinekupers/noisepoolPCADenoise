function [trl, Events] = ssmeg_ft_trial_fun(cfg, onsets, conditions)

% read the header information and the events from the data
hdr   = ft_read_header(cfg.dataset);

% Get onset of trigger events
trigger = find(onsets);

% determine the number of samples before and after the trigger
pretrig  = -round(cfg.trialdef.pre  * hdr.Fs);
posttrig =  round(cfg.trialdef.post * hdr.Fs);

for j = 1:length(trigger);  
      
    trlbegin = trigger(j) +  pretrig;     
    trlend   = trigger(j) + posttrig;       
    trl(j,:) = [trlbegin trlend pretrig];
    
    Events(j).trigger = trigger(j);
%     Events(j).channel = trigger(j).type; % I don't think we need this..
    Events(j).value   = conditions(j);
    Events(j).sample = 'trial';
    Events(j).type = 'trial';
end


return