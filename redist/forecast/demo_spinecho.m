% Spin echo examples for the forecast simulation framework
% 1. Air bubble
% 2. Titanium cylinder
% 3. Water-fat shift
%
% Note that these are simply literal copies of the gradient echo examples
% with the acquisition method changed to spin echo. The spatial deformation
% will be the same as in the gradient echo examples, while these is no
% dephasing.
%
% For spin echo simulation an additional k-space sampling time matrix is
% required: acquisition.kspaceSamplingTimesRefocused
% These sampling times are applied only to the delta-B0 term of the signal
% formula (which refocuses), while the original sampling times are applied
% only to the T2 term (which does not refocus).
% In most spin echo acquisitions the calculation will be:
% kspaceSamplingTimesRefocused = kspaceSamplingTimes - echoTime
% (i.e. the range of values will be -samplingInterval/2 to
% samplingInterval/2)

clear all
close all

%% 3D simulation of air bubble: 128x128x128 model, 64x64x64 acquisition k-space

Nmodel = 128;
Nacq = 64;
readoutDirection = 'z';
echoTime = 4;
samplingInterval = 8;

fprintf('========\n');
fprintf('Example 1: Air bubble in 3D\nModel matrix: %dx%dx%d\nAcquisition matrix: %dx%dx%d\nReadout direction: %s\nTE: %f ms\nSampling interval: %f ms\n', Nmodel, Nmodel, Nmodel, Nacq, Nacq, Nacq, readoutDirection, echoTime, samplingInterval);


% Create a spherical mask for air
[x,y,z] = meshgrid(linspace(-1,1,Nmodel),linspace(-1,1,Nmodel),linspace(-1,1,Nmodel));
airMask = sqrt(x.^2 + y.^2 + z.^2) <= 0.2;

% Create simulation model: Density 1 in background, density 0 inside air bubble
model = struct();
model.protonDensity = ones(Nmodel,Nmodel,Nmodel);
model.protonDensity(airMask) = 0;

susWater = -9e-6;
susAir = 0;

susceptibility = susWater * ones(Nmodel,Nmodel,Nmodel);
susceptibility(airMask) = susAir;

% Calculate susceptibility relative to background susceptibility to avoid boundary artifacts
susceptibility = susceptibility - susWater;

model.resolution = [Nacq/Nmodel Nacq/Nmodel Nacq/Nmodel];
fprintf('Calculating susceptibility-induced field shift... ');
tic
model.deltaB0 = calculateFieldShift(susceptibility, model.resolution);
toc

% Acquisition parameters
acquisition = struct();
[acquisition.kspaceSamplingTimes, acquisition.kspaceSamplingTimesRefocused] = calculateCartesianSamplingTimes('spinecho', [Nacq Nacq Nacq], readoutDirection, echoTime, samplingInterval);
acquisition.resolution = [1 1 1];

% Simulation
fprintf('Simulating... ');
tic
kspace = forecast(model, acquisition, struct('verbose', true));
toc
image = ifftc(kspace);


figure
subplot(2,2,1), imagesc(rot90(squeeze(model.protonDensity(:,fftCenter(Nmodel),:)),-1))
colormap(gray(256))
axis equal
title('3D Air bubble - Proton density')
subplot(2,2,2), imagesc(rot90(squeeze(model.deltaB0(:,fftCenter(Nmodel),:)),-1))
colormap(gray(256))
axis equal
title('Delta-B0')

subplot(2,2,3), imagesc(rot90(abs(squeeze(image(:,fftCenter(Nacq),:))),-1),[0 1.5])
colormap(gray(256))
axis equal
title('Z readout - Coronal slice - Magnitude')
subplot(2,2,4), imagesc(rot90(angle(squeeze(image(:,fftCenter(Nacq),:))),-1),[-pi pi])
colormap(gray(256))
axis equal
title('Phase')
drawnow


%% 3D simulation of titanium cylinder: 256x256x256 model, 128x128x128 acquisition k-space

Nmodel = 256;
Nacq = 128;
readoutDirection = 'z';
echoTime = 3;
samplingInterval = 2;

fprintf('========\n');
fprintf('Example 2: Titanium cylinder in 3D\nModel matrix: %dx%dx%d\nAcquisition matrix: %dx%dx%d\nReadout direction: %s\nTE: %f ms\nSampling interval: %f ms\n', Nmodel, Nmodel, Nmodel, Nacq, Nacq, Nacq, readoutDirection, echoTime, samplingInterval);


% Create a cylindrical mask
[x,y,z] = meshgrid(linspace(-1,1,Nmodel),linspace(-1,1,Nmodel),linspace(-1,1,Nmodel));
cylinderMask = sqrt(x.^2 + y.^2) < 0.1 & abs(z) <= 0.3;

% Create simulation model: Density 1 in background, density 0 inside cylinder
model = struct();
model.protonDensity = ones(Nmodel,Nmodel,Nmodel);
model.protonDensity(cylinderMask) = 0;

susWater = -9e-6;
susTitanium = 180e-6;

susceptibility = susWater * ones(Nmodel,Nmodel,Nmodel);
susceptibility(cylinderMask) = susTitanium;

% Calculate susceptibility relative to background susceptibility to avoid boundary artifacts
susceptibility = susceptibility - susWater;

model.resolution = [Nacq/Nmodel Nacq/Nmodel Nacq/Nmodel];
fprintf('Calculating susceptibility-induced field shift... ');
tic
model.deltaB0 = calculateFieldShift(susceptibility, model.resolution);
toc

% Acquisition parameters
acquisition = struct();
[acquisition.kspaceSamplingTimes, acquisition.kspaceSamplingTimesRefocused] = calculateCartesianSamplingTimes('spinecho', [Nacq Nacq Nacq], readoutDirection, echoTime, samplingInterval);
acquisition.resolution = [1 1 1];

% Simulation
fprintf('Simulating (takes ~1 minute)... ');
tic
kspace = forecast(model, acquisition);
toc
image = ifftc(kspace);


figure
subplot(2,2,1), imagesc(rot90(squeeze(model.protonDensity(:,fftCenter(Nmodel),:)),-1))
colormap(gray(256))
axis equal
title('3D Titanium cylinder - Proton density')
subplot(2,2,2), imagesc(rot90(squeeze(model.deltaB0(:,fftCenter(Nmodel),:)),-1))
colormap(gray(256))
axis equal
title('Delta-B0')

subplot(2,2,3), imagesc(rot90(abs(squeeze(image(:,fftCenter(Nacq),:))),-1),[0 1.5])
colormap(gray(256))
axis equal
title('Z readout - Coronal slice - Magnitude')
subplot(2,2,4), imagesc(rot90(angle(squeeze(image(:,fftCenter(Nacq),:))),-1),[-pi pi])
colormap(gray(256))
axis equal
title('Phase')
drawnow


%% 2D simulation of basic water-fat shift: 256x256 model, 128x128 acquisition k-space
% Note that this example assumes that a voxel is either completely water or
% completely fat. If you want to simulate a model with partial water/fat
% volumnes, define a separate water and fat matrix, do two separate
% simulations and add up the resulting k-spaces.
% Alse note that no susceptibility effects are calculated here (i.e. the
% model represents a 2D slice through an infinite cylinder along B0).

Nmodel = 256;
Nacq = 128;
echoTime = 4.4; % Note that changing the echo time to 2.2 ms does not create an out of phase effect, contrary to gradient echo
samplingInterval = 8;

fprintf('========\n');
fprintf('Example 3: Water-fat shift in 2D\nModel matrix: %dx%d\nAcquisition matrix: %dx%d\nTE: %f ms\nSampling interval: %f ms\n', Nmodel, Nmodel, Nacq, Nacq, echoTime, samplingInterval);


% Create a circular mask for fat
[x,y] = meshgrid(linspace(-1,1,Nmodel),linspace(-1,1,Nmodel));
waterMask = sqrt((x).^2 + y.^2) > 0.2;
fatMask = sqrt(x.^2 + (y).^2) <= 0.2;

% Create model: density 0.5 for water, density 1 for fat
model = struct();
model.resolution = [Nacq/Nmodel Nacq/Nmodel 1];

model.protonDensity = zeros(Nmodel,Nmodel);
model.protonDensity(waterMask) = 0.5;
model.protonDensity(fatMask) = 1;

% Convert fat off-resonance (220 Hz at 1.5T) to a B0 offset
gm = 42.57747892e6; % Hz / T
model.deltaB0 = zeros(Nmodel,Nmodel);
model.deltaB0(fatMask) = 220 / gm;

% Acquisition parameters - Readout in X dimension
acquisition = struct();
[acquisition.kspaceSamplingTimes, acquisition.kspaceSamplingTimesRefocused] = calculateCartesianSamplingTimes('spinecho', [Nacq Nacq 1], 'x', echoTime, samplingInterval);
acquisition.resolution = [1 1 1];

% Simulation
fprintf('X readout: Simulating... ');
tic
kspace = forecast(model, acquisition);
toc

figure
subplot(3,2,1), imagesc(model.protonDensity,[0 1.5])
colormap(gray(256))
axis equal
title('Water-fat shift - Proton density')
subplot(3,2,2), imagesc(model.deltaB0)
colormap(gray(256))
axis equal
title('W-F shift converted to delta-B0')

subplot(3,2,3), imagesc(abs(ifftc(kspace)),[0 1.5])
colormap(gray(256))
axis equal
title('Readout X - Magnitude')
subplot(3,2,4), imagesc(angle(ifftc(kspace)),[-pi pi])
colormap(gray(256))
axis equal
title('Phase')


% Acquisition parameters - Readout in Y dimension
acquisition = struct();
[acquisition.kspaceSamplingTimes, acquisition.kspaceSamplingTimesRefocused] = calculateCartesianSamplingTimes('spinecho', [Nacq Nacq 1], 'y', echoTime, samplingInterval);
acquisition.resolution = [1 1 1];

% Simulation
fprintf('Y readout: Simulating... ');
tic
kspace = forecast(model, acquisition);
toc

subplot(3,2,5), imagesc(abs(ifftc(kspace)),[0 1.5])
colormap(gray(256))
axis equal
title('Readout Y - Magnitude')
subplot(3,2,6), imagesc(angle(ifftc(kspace)),[-pi pi])
colormap(gray(256))
axis equal
title('Phase')
drawnow