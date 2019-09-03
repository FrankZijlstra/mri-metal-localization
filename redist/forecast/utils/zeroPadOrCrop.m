function [out] = zeroPadOrCrop (im, sz, padValue)
% Either zero-pads or crops matrix im to size sz, depending on whether the
% size increases of decreases

if (nargin < 3)
    padValue = 0;
end

if (numel(sz) == 2)
    sz(3) = 1;
end

out = zeroPad(crop(im, min(sz, [size(im,1) size(im,2) size(im,3)])), sz, padValue);

end
