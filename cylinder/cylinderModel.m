function [PD, dB, T2] = cylinderModel (xSpace, ySpace, zSpace, R, voxelSize, B0)

[x,y,z] = meshgrid(xSpace,ySpace,zSpace);

susTitanium = 182e-6; %https://www-mrsrl.stanford.edu/studygroup/2/Files/Schenck_susceptibility.pdf
susBackground = -9e-6;

if (~isempty(R))
    ps = [x(:) y(:) z(:)] * R;
end

PD = ones(size(x));
susceptibility = ones(size(x))*susBackground;

if (~isempty(R))
    % Smooth heaviside: 1/(1 + exp(-2*k*x))
    k = 20;

    % Cylinder
    dist1 = sqrt(ps(:,2).^2 + ps(:,1).^2) - 25/2;
    dist2 = abs(ps(:,3)) - 85/2;

    dist = max(dist1,dist2);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*H + (1-H).*susTitanium;
    PD = PD.*H;
end

% Calculate field-shift
dB = B0 * calculateFieldShift(susceptibility - susBackground, voxelSize);
T2 = Inf;

% Suppress spins outside of excitation bandwidth
gm = 42.57747892e6; % Hz / T
PD(abs(dB*gm)>6800/2) = 0;

end
