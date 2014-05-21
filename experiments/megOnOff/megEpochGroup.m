function epochGroup = megEpochGroup(okEpochs,group_epoch,shift_epoch)
% create epochGroup for doing denoising 

if notDefined('group_epoch'),   group_epoch  = 1;  end
if notDefined('shift_epoch'),   shift_epoch = 0; end

fprintf('(megEpochGroup) grouping everying %d epochs for denoising\n', group_epoch);
epochGroup = repmat(1:length(okEpochs)/group_epoch,group_epoch,1);
epochGroup = epochGroup(:);
epochGroup = circshift(epochGroup,shift_epoch);
epochGroup(1 : shift_epoch) = 0;
epochGroup(end-shift_epoch+1 : end) = 0;
epochGroup = epochGroup(okEpochs);
