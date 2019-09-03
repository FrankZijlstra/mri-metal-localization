function [data, position, orientation] = detectSeed(data)
% Detect next best seed, position and orientation are left empty if nothing suitable can be found

maskIndices = find(data.mask);
if (isempty(maskIndices))
    position = [];
    orientation = [];
    return;
end

%% Sort by regression error, correct the phase of the image, and repeat regression
[~,index] = sort(data.regressionError(maskIndices), 'ascend');

index = maskIndices(index(1));
[y,x,z] = ind2sub(size(data.image), index);
position = data.worldMatrix(1:3,:) * [x-1; y-1; z-1; 1];
orientation = data.library.templates(data.templateId(index)).orientationVector.';

data.image = correctPhase(data.image, index, data.library.templates(data.templateId(index)).image);
data.mask(index) = 0;
data.mask(getNeighbourhoodIndices(size(data.image), index, ones(3,3,3))) = 0;

re = matchLibraryRegression(data.image, data.templateId, data.library, data.mask);
data.regressionError(data.mask) = re(data.mask);

end
