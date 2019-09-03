function [out] = crop (im, sz)
% Crops matrix im to size sz

if (length(sz) == 2)
    sz(3) = 1;
end

startX = fftCenter(size(im,1)) - fftCenter(sz(1)) + 1;
startY = fftCenter(size(im,2)) - fftCenter(sz(2)) + 1;
startZ = fftCenter(size(im,3)) - fftCenter(sz(3)) + 1;

out = im(startX:startX+sz(1)-1,startY:startY+sz(2)-1,startZ:startZ+sz(3)-1);

end
