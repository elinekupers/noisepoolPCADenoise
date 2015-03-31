function chanNum0 = megGetOrigChannel(chanNum,badChannels,orig2new)
% mapping between channel indices between original channel space and the
% new space with bad channels discarded
if notDefined('orig2new'), orig2new = true; end

if orig2new %figure out index in the vector with badChannels discarded
    chanNum0 = nan(size(chanNum));
    for jj = 1:length(chanNum)
        tmp = zeros(1,157);
        tmp(chanNum(jj))=1;
        if ~isempty(find(tmp(~badChannels), 1))
            chanNum0(jj)= find(tmp(~badChannels));
        end
    end
else % figure out index in original space
    tmp = zeros(1,157);
    tmp(~badChannels) = 1:sum(~badChannels);
    chanNum0 = zeros(size(chanNum));
    for jj = 1:length(chanNum)
        chanNum0(jj) = find(tmp==chanNum(jj));
    end
end
