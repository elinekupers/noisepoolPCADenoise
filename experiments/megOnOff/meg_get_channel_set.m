function channel_nums = meg_get_channel_set(channel_group)
%
% get a list of MEG channel numbers that belong to a group of type
% 'channel_type', for example 'back_half' which are all channels in the
% back half of the head
%
% Inputs
%   channel_type: string indicating the grouping, e.g., 'back_half'
%
% Ouputs
%   channel_nums: vector of channel numbers belonging to the group
%
% Example
%   channel_nums = meg_get_channel_set('back_half');

if isempty(which('hdr.mat')), error('denoiseproject needed on path'); end

% hdr is a standard header file for our MEG 160 system. It contains basic
% information about spatial configuration of the channels, sampling
% frequency, and so forth. in principle, there should be a hdr file for
% every experiment run, and only the appropriate hdr file should be loaded.
% in this case, we use a default file because we want a simple piece of
% information that should not change with experiment (though it will change
% if a different MEG system is used).
load('hdr')

cfg = [];
cfg.interpolation = 'nearest';
layout = ft_prepare_layout(cfg, hdr);
num_channels = length(layout.cfg.channel);

switch lower(channel_group)
    case {'back_half' 'back'}
        channel_nums = find(layout.pos(:,2)<0);
    case {'front_half' 'front'}
        channel_nums = find(layout.pos(:,2)>0);
    otherwise 
        error('unknown channel type %s', channel_group)
end

channel_nums = channel_nums(channel_nums <= num_channels);
return