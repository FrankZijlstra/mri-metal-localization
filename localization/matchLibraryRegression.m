function [regressionError, regressionSlope, regressionOffset] = matchLibraryRegression (image, templateId, library, mask)
% After initial matching with matchLibrary, this can be used to calculate
% a local least squares matching error of a linear regression of the signal
% in the image with the template signal.
%
% The regression is only applied within the specified mask. If mask is
% unspecified, it is applied to all voxels (this is very slow!).

differenceThreshold = 0.2; % Match only part of the template that changed at least this percentage from the background intensity

if (nargin < 4)
    mask = true(size(image)); % SLOW!
end

inds = find(mask);

templateInds = templateId(inds);
uniqueTemplateInds = unique(templateInds);

regressionError = zeros(size(image));
regressionSlope = zeros(size(image));
regressionOffset = zeros(size(image));

for I=1:length(uniqueTemplateInds)
    IND = uniqueTemplateInds(I);
    
    templateImage = library.templates(IND).image;
    regMask = abs(templateImage - library.background)>differenceThreshold*max(abs(library.background(:))); % Threshold is relative to max magnitude in background

    tmpInds = find(templateInds == IND);

    for J=1:length(tmpInds)
        [tx,ty,tz] = ind2sub(size(image),inds(tmpInds(J)));
        
        startX = tx-fftCenter(size(templateImage,1))+1;
        startY = ty-fftCenter(size(templateImage,2))+1;
        startZ = tz-fftCenter(size(templateImage,3))+1;

        if(startX <= 0 || startY <= 0 || startZ <= 0 || startX+size(templateImage,1)-1>size(image,1) || startY+size(templateImage,2)-1>size(image,2) || startZ+size(templateImage,3)-1>size(image,3))
            mask(inds(tmpInds(J))) = false;
            continue;
        end

        scanImage = image(startX:startX+size(templateImage,1)-1,startY:startY+size(templateImage,2)-1,startZ:startZ+size(templateImage,3)-1);

        % Find linear mapping from scanImage to templateImage and calculate
        % residual error        
        R = [scanImage(regMask) ones(size(scanImage(regMask)))] \ templateImage(regMask);
        err = sum(abs(R(1) * scanImage(regMask) + R(2) - templateImage(regMask)).^2);

        regressionError(tx,ty,tz) = err;
        regressionSlope(tx,ty,tz) = R(1);
        regressionOffset(tx,ty,tz) = R(2);
    end
end

regressionError(~mask) = max(regressionError(:));

end

