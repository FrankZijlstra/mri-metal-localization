function [PD, dB, T2] = seedModel_Best2301 (xSpace, ySpace, zSpace, R, voxelSize, B0)

[x,y,z] = meshgrid(xSpace,ySpace,zSpace);

susTitanium = 1.795e-4; % Wachowicz
susTungsten = 59e-6;
susAir = 0.36e-6;
susPolystyrene = -8.21e-6; % http://pubs.rsc.org/en/content/articlehtml/2012/an/c2an35199d
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
    dist = sqrt((2*(ps(:,3)-2.2750)).^2 + ps(:,2).^2 + ps(:,1).^2) - 0.4 + 1000*(abs(ps(:,3))<=2.2750);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susTitanium;
    PD = PD.*H;

    dist = sqrt((2*(ps(:,3)+2.2750)).^2 + ps(:,2).^2 + ps(:,1).^2) - 0.36 + 1000*(abs(ps(:,3))<=2.2750);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susTitanium;
    PD = PD.*H;

    dist = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.4 + 1000*(abs(ps(:,3))>2.2750);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susTitanium;
    PD = PD.*H;


    dist = sqrt((2*(ps(:,3)-2.2750)).^2 + ps(:,2).^2 + ps(:,1).^2) - 0.36 + 1000*(abs(ps(:,3))<=2.2750);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susAir;
    PD = PD.*H;

    dist = sqrt((2*(ps(:,3)+2.2750)).^2 + ps(:,2).^2 + ps(:,1).^2) - 0.32 + 1000*(abs(ps(:,3))<=2.2750);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susAir;
    PD = PD.*H;

    dist = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.32 + 1000*(abs(ps(:,3))>2.2750);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susAir;
    PD = PD.*H;

    dist = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.225 + 1000*(abs(ps(:,3))>2.025);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susPolystyrene;
    PD = PD.*H;

    dist = sqrt(ps(:,2).^2 + ps(:,1).^2) - 0.125 + 1000*(abs(ps(:,3))>1.875);
    H = reshape(1./(1 + exp(-2*k*dist)), size(PD));
    susceptibility = susceptibility.*(H) + (1-H).*susTungsten;
    PD = PD.*H;
end

% Calculate field-shift
dB = B0 * calculateFieldShift(susceptibility - susBackground, voxelSize);
T2 = 50;

end
