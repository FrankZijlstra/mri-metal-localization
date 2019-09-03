function [ space ] = calculateSimulationCoordinates (fov, n)
% Calculate the coordinates of the voxels centers for a given field of view
% (given as the limits of the field of view) and number of voxels within
% that field of view. This function makes sure the center of the field of
% view is at the fft center (fftCenter(n)).
%
% Input:
%   fov: [min max]
%   n:   number of voxels within this FOV


%% Calculate the voxel centers of the first and last voxel.
% This corrects for both the voxel size and the possible asymmetry created
% by the fft center.
voxelSize = (fov(2) - fov(1))/n;

if (mod(n,2) == 0)
    fov(1) = fov(1);
    fov(2) = fov(2) - voxelSize;
else
    fov(1) = fov(1) + voxelSize/2;
    fov(2) = fov(2) - voxelSize/2;
end

%% Calculate coordinates
space = linspace(fov(1), fov(2), n);

end
