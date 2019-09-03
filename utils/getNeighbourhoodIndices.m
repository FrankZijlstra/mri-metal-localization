function [inds] = getNeighbourhoodIndices (sz, index, nhood)

szn = size(nhood);

[y,x,z] = ind2sub(sz, index);
[y2,x2,z2] = ind2sub(szn, find(nhood));
C = [y+y2-1-floor(szn(1)/2) x+x2-1-floor(szn(2)/2) z+z2-1-floor(szn(3)/2)];
valid = all(C >= 1 & C <= repmat(sz,size(C,1),1),2);
inds = sub2ind(sz, C(valid,1), C(valid,2), C(valid,3));


end

