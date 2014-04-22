function newArr = to157chan(currArr,inds,how)
% newArr should be nconds x N, where is N less than 157
% inds is 1x157 boolean where N entries are 1's
%
if isnumeric(how)
    newArr = how*ones(size(currArr,1),157);
elseif strcmp(how,'nans')
    newArr = nan(size(currArr,1),157);
elseif strcmp(how,'zeros')
    newArr = zeros(size(currArr,1),157);
else
    error('input "how" not recognized');
end
newArr(:,inds) = currArr;
