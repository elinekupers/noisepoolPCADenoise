function out = getsignalnoise(model,whichbeta,type)
% returns signal, noise, or snr from modelfit
% whichbeta: condition number, e.g., 1, or 1:3
% type: what to return, e.g., 'SNR'

if notDefined('type'), type = 'SNR'; end
if notDefined('whichbeta'), whichbeta = 1:size(model.beta,1); end

signal = max(abs(model.beta_md(whichbeta,:)),[],1);
noise  = mean(model.beta_se(whichbeta,:),1);
switch type
    case 'S'
        out = signal;
    case 'N'
        out = noise;
    case 'SNR'
        out = signal./noise;
end
