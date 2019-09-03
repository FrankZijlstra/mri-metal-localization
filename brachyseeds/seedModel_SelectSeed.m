function [PD, dB, T2] = seedModel_SelectSeed (xSpace, ySpace, zSpace, R, voxelSize, B0)

[x,y,z] = meshgrid(xSpace,ySpace,zSpace);

susTit = 1.795e-4; % Wachowicz
susAg = -2.38e-5; % Wachowicz
susAir = 0.36e-6;
susBackground = -9e-6;

if (~isempty(R))
    ps = [x(:) y(:) z(:)] * R;
end

PD = ones(size(x));
susceptibility = ones(size(x))*susBackground;

% Smooth heaviside: 1/(1 + exp(-2*k*x))
% k = 40 good for modelMultiplier = 16
k = 40;

if (~isempty(R))
    % Titanium
    distC1 = sqrt((ps(:,3)-1.8).^2 + ps(:,2).^2 + ps(:,1).^2) - 0.4;
    distC2 = sqrt((ps(:,3)+1.8).^2 + ps(:,2).^2 + ps(:,1).^2) - 0.4;
    distCy = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.4 + 1000*(abs(ps(:,3))>1.8);

    dist = min(min(distC1,distC2),distCy);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*H + (1-H).*susTit;
    PD = PD.*H;

    % Air
    dist = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.35 + 1000*(abs(ps(:,3))>1.7);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*H + (1-H).*susAir;
    PD = PD.*H;

    % Silver
    dist = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.25 + 1000*(abs(ps(:,3))>1.6);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*H + (1-H).*susAg;
    PD = PD.*H;
end

% Calculate field-shift
dB = B0 * calculateFieldShift(susceptibility - susBackground, voxelSize);
T2 = 50;


end
