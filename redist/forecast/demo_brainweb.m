% Simulation of a brainweb dataset
% This script will automatically download and unzip a binary brain
% segmentation from brainweb (http://www.bic.mni.mcgill.ca/brainweb/).
% Note that simulation takes about 5 minutes on a modern intel CPU and
% requires around 6.5G of RAM.
%
% Included effects:
%   T1 (Analytically)
%   T2
%   Susceptibility
%   Water-fat shift (voxels are either 100% water or 100% fat)

clear all
close all

%% Parameters

echoTime = 15;
repetitionTime = 20;
samplingInterval = 4;
flipAngle = 30;

readoutDirection = 'x';

B0 = 1.5; % Tesla (note that tissue parameters below are hardcoded for 1.5T)


%% Read raw brainweb data

% If file doesn't exist, download from brainweb: http://www.bic.mni.mcgill.ca/brainweb/
if (~exist('subject04_crisp_v.rawb', 'file'))
    if (verLessThan('matlab','8.5'))
        error('Matlab R2014B and earlier cannot download brainweb dataset automatically. Please download the dataset from http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?alias=subject04_crisp&download=1 in "raw byte" format and place the unpacked subject04_crisp_v.rawb file in the current directory.');
    end
    
    websave('subject04_crisp_v.rawb.gz', 'http://brainweb.bic.mni.mcgill.ca/cgi/brainweb1?alias=subject04_crisp&download=1', 'do_download_alias', 'subject04_crisp', 'format_value', 'raw_byte', 'zip_value', 'gnuzip', weboptions('RequestMethod','post'));
    gunzip('subject04_crisp_v.rawb.gz');
end

if (~exist('subject04_crisp_v.rawb', 'file'))
    error('Failed to download brainweb dataset?');
end

% Read raw bytes
fp = fopen('subject04_crisp_v.rawb', 'rb');
data = fread(fp,'uint8=>uint8');
fclose(fp);

% Reshape and flip dimensions around a bit
data = reshape(data,[362 434 362]);
data = flip(permute(data,[2 1 3]),1);


%% Make model out of brainweb data
% Note that we guessed some of the tissue parameters. It doesn't make much
% difference for the simulation, but if you have better values feel free to
% let us know!

susWater = -9.04e-6;
susAir = 0.40e-6;
susGM = -8.97e-6;
susWM = -8.80e-6;
susCSF = susWater; % guessed
susFat = -7.79e-6;
susMuscle = susWater; % guessed
susSkin = susWater; % guessed
susBone = -8.44e-6;
susBlood = -9.12e-6; % Oxygenated blood
susConnectiveTissue = susWater; % guessed
susDura = susWater; % guessed
susMarrow = susBone; % guessed

susceptibility = zeros(size(data));
susceptibility(data == 0) = susAir;
susceptibility(data == 1) = susCSF;
susceptibility(data == 2) = susGM;
susceptibility(data == 3) = susWM;
susceptibility(data == 4) = susFat;
susceptibility(data == 5) = susMuscle;
susceptibility(data == 6) = susSkin;
susceptibility(data == 7) = susBone;
susceptibility(data == 8) = susBlood;
susceptibility(data == 9) = susConnectiveTissue;
susceptibility(data == 10) = susDura;
susceptibility(data == 11) = susMarrow;

protonDensity = zeros(size(data));
protonDensity(data == 0) = 0; % Air
protonDensity(data == 1) = 1; % CSF
protonDensity(data == 2) = 0.8; % GM
protonDensity(data == 3) = 0.65; % WM
protonDensity(data == 4) = 0.9; % Fat
protonDensity(data == 5) = 0.7; % Muscle (guessed)
protonDensity(data == 6) = 0.7; % Skin (guessed)
protonDensity(data == 7) = 0.9; % Skull (guessed)
protonDensity(data == 8) = 1; % Vessels
protonDensity(data == 9) = 0.7; % Connective tissue (guessed)
protonDensity(data == 10) = 0.7; % Dura (guessed)
protonDensity(data == 11) = 0.8; % Marrow (guessed)

T1 = zeros(size(data));
T1(data == 0) = 0; % Air
T1(data == 1) = 4500; % CSF
T1(data == 2) = 950; % GM
T1(data == 3) = 600; % WM
T1(data == 4) = 250; % Fat
T1(data == 5) = 900; % Muscle
T1(data == 6) = 900; % Skin (guessed)
T1(data == 7) = 500; % Skull (guessed)
T1(data == 8) = 1200; % Vessels
T1(data == 9) = 900; % Connective tissue (guessed)
T1(data == 10) = 900; % Dura (guessed)
T1(data == 11) = 500; % Marrow (guessed)

T2 = zeros(size(data));
T2(data == 0) = 0; % Air
T2(data == 1) = 2200; % CSF
T2(data == 2) = 100; % GM
T2(data == 3) = 80; % WM
T2(data == 4) = 60; % Fat
T2(data == 5) = 50; % Muscle
T2(data == 6) = 50; % Skin (guessed)
T2(data == 7) = 1; % Skull (guessed)
T2(data == 8) = 200; % Vessels
T2(data == 9) = 50; % Connective tissue (guessed)
T2(data == 10) = 50; % Dura (guessed)
T2(data == 11) = 8; % Marrow (guessed)


% Model structure
model = struct();
model.resolution = [0.5 0.5 0.5];
model.T2 = T2;

% Calculate the steady-state analytically
E1 = exp(-repetitionTime ./ T1);
model.protonDensity = protonDensity .* sin(flipAngle * pi/180) .* (1 - E1) ./ (1 - E1 .* cos(flipAngle * pi/180));

% Calculate delta-B0
fprintf('Calculating susceptibility-induced field shift... ');
tic
model.deltaB0 = B0 * calculateFieldShift(susceptibility - susAir, model.resolution);
toc

gm = 42.57747892e6; % Hz / T

% Add fat off-resonance at 1.5T, converted from Hz to Tesla
model.deltaB0(data == 4) = model.deltaB0(data == 4) + 220/gm;


% Acquisition structure
acquisition = struct();
acqSizeYXZ = size(model.protonDensity)/2;
acquisition.kspaceSamplingTimes = calculateCartesianSamplingTimes('gradientecho', acqSizeYXZ([2 1 3]), readoutDirection, echoTime, samplingInterval);
acquisition.resolution = [1 1 1];

%% Simulation
fprintf('Simulating... ');
tic
kspace = forecast(model, acquisition, struct('verbose', true));
toc
image = ifftc(kspace);

%% Display results
figure
subplot(2,3,1), imagesc(abs(image(:,:,fftCenter(size(image,3)))),[0 0.125])
colormap(gray(256))
axis equal
subplot(2,3,2), imagesc(rot90(abs(squeeze(image(:,fftCenter(size(image,2))+5,:)))),[0 0.125])
colormap(gray(256))
axis equal
subplot(2,3,3), imagesc(rot90(abs(squeeze(image(fftCenter(size(image,1)),:,:)))),[0 0.125])
colormap(gray(256))
axis equal
subplot(2,3,4), imagesc(angle(image(:,:,fftCenter(size(image,3)))),[-pi pi])
colormap(gray(256))
axis equal
subplot(2,3,5), imagesc(rot90(angle(squeeze(image(:,fftCenter(size(image,2))+5,:)))),[-pi pi])
colormap(gray(256))
axis equal
subplot(2,3,6), imagesc(rot90(angle(squeeze(image(fftCenter(size(image,1)),:,:)))),[-pi pi])
colormap(gray(256))
axis equal

