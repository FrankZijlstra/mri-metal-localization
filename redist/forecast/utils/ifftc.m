function res = ifftc(x, siz)
% N-D Inverse centered FFT

if (nargin < 2)
    res = sqrt(numel(x))*fftshift(ifftn(ifftshift(x)));
else
    res = sqrt(prod(siz))*fftshift(ifftn(ifftshift(x),siz));
end

end
