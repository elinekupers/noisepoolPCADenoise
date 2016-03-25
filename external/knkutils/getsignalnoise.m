function out = getsignalnoise(model,whichbeta,type)
% returns signal, noise, or snr from modelfit
% whichbeta: condition number, e.g., 1, or 1:3
% type: what to return, e.g., 'SNR'

if notDefined('type'), type = 'SNR'; end
if notDefined('whichbeta'), whichbeta = 1:size(model.beta,1); end

% Combine boots of channels if using Elekta system
if size(model.beta,2) > 200
    beta = mean([model.beta(whichbeta,1:2:end,:);model.beta(whichbeta,2:2:end,:)],1);
     
    temp = prctile(beta,[16 50 84],3);
    model.beta_md = temp(:,:,2);
    model.beta_se = diff(temp(:,:,[1 3]),[],3)/2;
    
    signal = max(model.beta_md,[],1);
    noise  = mean(model.beta_se,1);
else


    % signal = max(abs(model.beta_md(whichbeta,:)),[],1); % old way of
    % computing Signal
    
    % No absolute value of signal
    signal = max(model.beta_md(whichbeta,:),[],1);
    noise  = mean(model.beta_se(whichbeta,:),1);
end


switch type
    case 'S'
        out = signal;
    case 'N'
        out = noise;
    case 'SNR'
        out = signal./noise;
end
