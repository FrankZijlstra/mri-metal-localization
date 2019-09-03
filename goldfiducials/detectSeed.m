function [data, position, orientation] = detectSeed(data)
% Detect next best seed, position and orientation are left empty if nothing suitable can be found

maskIndices = find(data.mask);
if (isempty(maskIndices))
    position = [];
    orientation = [];
    return;
end

%% Sort by regression error and only accept detections that improve phase smoothness after correction, afterwards, correct the phase of the image, and repeat regression
[~,index] = sort(data.regressionError(maskIndices), 'ascend');

found = false;
for I=1:length(index)
    smooth = correctedPhaseSmoothness(data.image, maskIndices(index(I)), data.library.templates(data.templateId(maskIndices(index(I)))).image);

    if (smooth > 0)
        found = true;
        break;
    end
end

if (~found)
    position = [];
    orientation = [];
    return;
end

index = maskIndices(index(I));
[y,x,z] = ind2sub(size(data.image), index);
position = data.worldMatrix(1:3,:) * [x-1; y-1; z-1; 1];
orientation = data.library.templates(data.templateId(index)).orientationVector.';

data.image = correctPhase(data.image, index, data.library.templates(data.templateId(index)).image);
data.mask(getNeighbourhoodIndices(size(data.image), index, ones(3,3,3))) = 0;

re = matchLibraryRegression(data.image, data.templateId, data.library, data.mask);
data.regressionError(data.mask) = re(data.mask);


end
