function newArray = to157chan(currentArray,inds,padval)
% currentArray should be nconds x N, where N is less than 157
% inds is 1x157 boolean where N entries are 1's
%
if isnumeric(padval)
    newArray = padval*ones(size(currentArray,1),157);
elseif strcmp(padval,'nans')
    newArray = nan(size(currentArray,1),157);
elseif strcmp(padval,'zeros')
    newArray = zeros(size(currentArray,1),157);
else
    error('input "padval" not recognized');
end
newArray(:,inds) = currentArray;
