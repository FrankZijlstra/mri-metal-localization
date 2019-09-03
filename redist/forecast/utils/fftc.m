function res = fftc(x, siz)
% N-D centered FFT

if (nargin < 2)
    res = 1/sqrt(numel(x))*fftshift(fftn(ifftshift(x)));
else
    res = 1/sqrt(prod(siz))*fftshift(fftn(ifftshift(x), siz));
end

end