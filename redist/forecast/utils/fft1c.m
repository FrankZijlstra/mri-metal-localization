function res = fft1c(x, n, dim)
% 1D centered FFT
% TODO: parameter n has unexpected results due to fft shifts

if (nargin < 3 || isempty(dim))
    dim = find(size(x)~=1,1);
end

if (nargin < 2 || isempty(n))
    n = size(x,dim);
end

res = 1/sqrt(n)*fftshift(fft(ifftshift(x,dim),n,dim),dim);

end
