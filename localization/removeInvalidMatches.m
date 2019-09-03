function [match] = removeInvalidMatches(match, library)
% Mask out matches where the template overlapped with the edge of the image

sz = size(library.background);
fc = fftCenter(sz);
beginP = sz - fc + 1 + 1; % One voxel less than actual border, no idea why
endP = fc - 1 + 1; % One voxel less than actual border, no idea why

valid = zeros(size(match));
valid(beginP(1):end-endP(1),beginP(2):end-endP(2),beginP(3):end-endP(3)) = 1;

match = match.*valid;

end

