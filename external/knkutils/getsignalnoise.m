function out = getsignalnoise(model,whichContrasts,whichOutput,badChannels)
% returns signal, noise, or snr from modelfit
% whichbeta: condition number, e.g., 1, or 1:3
% type: what to return, e.g., 'SNR'
num_conditions = size(model.beta,1);
if notDefined('whichOutput'), whichOutput = 'SNR'; end
if notDefined('whichContrasts'), whichContrasts = eye(num_conditions); end
if notDefined('badChannels') 
    if size(model.beta,2) > 157; badChannels = zeros(1, 204);
    else, badChannels = zeros(1, 157); end; 
end

num_contrasts  = size(whichContrasts,1);
assert(isequal(size(whichContrasts, 2), num_conditions))

% % Combine boots of channels if using Elekta system
% if size(model.beta,2) > 157
%     if sum(badChannels) >= 1
%         % Identify bad channels in data
%         dataIn = nan(num_conditions,204,100);
%         dataIn(:,~badChannels,:) = model.beta(:,:,:); 
%         beta = nanmean([dataIn(:,1:2:end,:);dataIn(:,2:2:end,:)],1);
%     else
%         beta = mean([model.beta(:,1:2:end,:);model.beta(:,2:2:end,:)],1);
%     end
% else
%     beta = model.beta;
% 
% end

beta = model.beta;

% Signal is max of betas across conditions. Noise is mean of se across
% conditions.
sz = size(beta);
beta = reshape(beta, sz(1), sz(2)*sz(3));
beta = whichContrasts*beta;
beta = reshape(beta, num_contrasts, sz(2), sz(3));
temp = prctile(beta,[16 50 84],3);

model.beta_md = temp(:,:,2);
model.beta_se = diff(temp(:,:,[1 3]),[],3)/2;

signal = max(model.beta_md,[],1);
noise  = mean(model.beta_se,1);


switch whichOutput
    case 'S'
        out = signal;
    case 'N'
        out = noise;
    case 'SNR'
        out = signal./noise;
end
