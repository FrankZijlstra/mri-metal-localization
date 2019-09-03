function [match, templateId] = matchLibrary (image, library, kSpaceMask)
% Run Phase Correlation matching on the specified image using the templates
% from the specified library.
%
% image should match the orientation of the image of that was used to
% create the library

if (nargin < 3 || isempty(kSpaceMask))
    kSpaceMask = ones(size(image));
end

if (~all(size(kSpaceMask) == size(image)))
    error('Mask size must match image size');
end

match = zeros(size(image));
templateId = ones(size(image)) * -1;

F = fftc(image);

for I=1:length(library.templates)
    fprintf('Matching: %d\n', I);
    
    T = fftc(zeroPad(library.templates(I).image - library.background, size(image)));

    % Phase Correlation template matching
    M = F.*conj(T);
    M = M ./ (abs(M)+eps) .* kSpaceMask; % Condition could be changed to require F and T to both have sufficient intensity (i.e. reliable k-space phase information)
    
    M = ifftc(M);
    
    templateId(abs(M) > abs(match)) = I;
    match(abs(M) > abs(match)) = M(abs(M) > abs(match));
end

end
