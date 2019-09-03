clear all
close all

dicomDir = './testdata/brachyseeds/';
seriesDescription = 'mFFE';
echoNumber = 2;

% dicomDir = './testdata/goldfiducials/';
% seriesDescription = 't mFFE 3D Seeds';
% echoNumber = 2;

% dicomDir = '../testdata/cylinder/';
% seriesDescription = 'FFE 1mm 0 deg wfs1';
% echoNumber = 2;


% This test simulates an image with the letters AP/RL/FH on the
% corresponding axes in the image. This tests whether all coordinate
% conversions are working as intended.

%% Add necessary paths
addpath('./simulation')
addpath('./localization')
addpath('./redist/dicomseries-matlab/dicomseries')
addpath('./redist/forecast')
addpath('./redist/forecast/utils')

%% Set library creation parameters
modelParameters.objectFunction = @testObjectFunction;
modelParameters.objectOrientations = [1 0 0]; % TODO: Does not include roll
modelParameters.modelSize = []; % voxels, scanner AP RL FH, keep empty to use scan size times modelMultiplier
modelParameters.modelMultiplier = 1;
modelParameters.modelFOV = []; % mm, scanner AP RL FH, keep empty to use scan FOV


% Simulate only center of scan region
% Set to region size (in voxels) in MPS (Measurement, Phase, Slice),
% or set to a single scalar n if the region is cubic (equal to [n n n]).
% Set to false to disable.
modelParameters.smallFOVsize = []; % voxels, MPS

% Post processing functions after simulation
options.postProcessImage = @(x) conj(x); % Convert from left-handed to right-handed spin


%% Load complex image from dicoms
dicomdict('set', 'dicom-dict-philips.txt');

fprintf('Scanning dicom directory...\n');

partitions = readDicomSeries(dicomDir, struct('verbose', true, 'recursive', true));

fprintf('Reading images...\n');
[dicomImage, dicomInfo] = readDicomSeriesImage(dicomDir, partitions, struct('SeriesDescription', seriesDescription, 'ImageType', 'ORIGINAL\PRIMARY\M_FFE\M\FFE', 'EchoNumber', echoNumber));
dicomImage = rescaleDicomImage(dicomImage, dicomInfo);

%% Create "library"
imageParameters = getImageParametersFromDicomPhilips(dicomInfo);
library = createLibrary(modelParameters, imageParameters, options);


%% Show an overlay of the simulated image with AP/RL/FH to check if this is consistent with the acquired image.
tmp = zeroPad(library.templates(1).image, size(dicomImage));
tmp = abs(tmp) / max(abs(tmp(:)));

im = dicomImage;
im = abs(im) / max(abs(im(:)));

im = im + tmp;

figure, imshow(im(:,:,round(size(dicomImage,3)/2)), [])
figure, imshow(squeeze(im(:,round(size(dicomImage,2)/2),:)), [])
figure, imshow(squeeze(im(round(size(dicomImage,1)/2),:,:)), [])


