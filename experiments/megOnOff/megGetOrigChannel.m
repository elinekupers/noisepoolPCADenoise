function chanNum0 = megGetOrigChannel(chanNum,badChannels)
% figure out index in the vector with badChannels discarded

tmp = zeros(1,157);
tmp(chanNum)=1;
chanNum0= find(tmp(~badChannels));

