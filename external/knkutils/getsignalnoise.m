function out = getsignalnoise(model,whichContrasts,whichOutput,varargin)
% returns signal, noise, or snr from modelfit
% whichbeta: condition number, e.g., 1, or 1:3
% type: what to return, e.g., 'SNR'
num_conditions = size(model.beta,1);
if notDefined('whichOutput'), whichOutput = 'SNR'; end
if notDefined('whichContrasts'), whichContrasts = eye(num_conditions); end


num_contrasts  = size(whichContrasts,1);
assert(isequal(size(whichContrasts, 2), num_conditions))

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
