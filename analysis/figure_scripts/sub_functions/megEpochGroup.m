function epochGroup = megEpochGroup(okEpochs,group_epoch,shift_epoch,remove_unequal_group)
% create epochGroup for doing denoising 

if notDefined('group_epoch'),   group_epoch  = 1;  end
if notDefined('shift_epoch'),   shift_epoch = 0; end
if notDefined('remove_unequal_group'),  remove_unequal_group = false; end

fprintf('(megEpochGroup) grouping every other %d epochs for denoising\n', group_epoch);
epochGroup = repmat(1:length(okEpochs)/group_epoch,group_epoch,1);
epochGroup = epochGroup(:);
epochGroup = circshift(epochGroup,shift_epoch);
epochGroup(1 : shift_epoch) = 0;
epochGroup(end-shift_epoch+1 : end) = 0;
epochGroup = epochGroup(okEpochs);

% Remove any skipped numbers
[~, ~, epochGroup] = unique(epochGroup);

% ensure row vector
epochGroup =  epochGroup(:)';

% ensure groups have equal numbers of epochs
if remove_unequal_group
    for ii = 1:max(epochGroup)
        c = epochGroup==ii;
        if sum(c) < group_epoch
            epochGroup(c) = 0;
        end
    end
end