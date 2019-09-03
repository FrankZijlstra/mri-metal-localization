function [image] = correctPhase(image, index, template)
% Rewinds the phase of the template at position indicated by index

[tx,ty,tz] = ind2sub(size(image),index);

startX = tx-fftCenter(size(template,1))+1;
startY = ty-fftCenter(size(template,2))+1;
startZ = tz-fftCenter(size(template,3))+1;

if (startX <= 0 || startY <= 0 || startZ <= 0 || startX+size(template,1)-1>size(image,1) || startY+size(template,2)-1>size(image,2) || startZ+size(template,3)-1>size(image,3))
    % TODO: Could do partial correction, but realistically this never happens
    return;
end

image(startX:startX+size(template,1)-1,startY:startY+size(template,2)-1,startZ:startZ+size(template,3)-1) = image(startX:startX+size(template,1)-1,startY:startY+size(template,2)-1,startZ:startZ+size(template,3)-1) .* exp(-1i * angle(template));

end

