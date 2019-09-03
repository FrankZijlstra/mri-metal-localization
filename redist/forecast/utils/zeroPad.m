function [out] = zeroPad (im, sz, value)
% Zero-pads matrix im to size sz, optionally with specified value

if (nargin == 3)
    out = value * ones(sz);
else
    out = zeros(sz);
end

if (numel(sz) == 2)
    sz(3) = 1;
end

startX = fftCenter(sz(1)) - fftCenter(size(im,1)) + 1;
startY = fftCenter(sz(2)) - fftCenter(size(im,2)) + 1;
startZ = fftCenter(sz(3)) - fftCenter(size(im,3)) + 1;

out(startX:startX+size(im,1)-1,startY:startY+size(im,2)-1,startZ:startZ+size(im,3)-1) = im;

end
