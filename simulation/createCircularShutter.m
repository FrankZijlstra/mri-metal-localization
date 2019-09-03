function [ kSpaceMask ] = createCircularShutter (sz, dicomInfo)
% Creates circular k-space mask in the phase encoding directions

% dicomInfo tells us which of the first two directions is the phase
% encoding

if (strcmp(getDicomAttribute(dicomInfo, 'InPlanePhaseEncodingDirection'), 'ROW'))
    [y, x] = meshgrid(linspace(-1,1,sz(3)),linspace(-1,1,sz(2)));
    dist = sqrt(x.^2 + y.^2);
    kSpaceMask = dist <= 1;
    kSpaceMask = repmat(permute(kSpaceMask,[3 1 2]), [sz(1) 1 1]);
else
    [y, x] = meshgrid(linspace(-1,1,sz(3)),linspace(-1,1,sz(1)));
    dist = sqrt(x.^2 + y.^2);
    kSpaceMask = dist <= 1;
    kSpaceMask = repmat(permute(kSpaceMask,[1 3 2]), [1 sz(2) 1]);
end


end
