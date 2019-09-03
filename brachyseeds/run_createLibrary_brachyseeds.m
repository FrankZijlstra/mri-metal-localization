clear all
close all

dicomDir = '../testdata/brachyseeds/';
seriesDescription = 'mFFE';
echoNumber = 2;

%% Add necessary paths
addpath('../simulation')
addpath('../localization')
addpath('../redist/dicomseries-matlab/dicomseries')
addpath('../redist/forecast')
addpath('../redist/forecast/utils')

%% Set library creation parameters

% Increase modelMultiplier below to increase the resolution of the
% simulation (also takes longer). Increase the argument to generateRotations
% further below to increase the number of rotations that are simulated
% (takes a lot longer...).

modelParameters.objectFunction = @seedModel_Best2301;
modelParameters.objectOrientations = generateRotations(1);
modelParameters.modelSize = []; % voxels, scanner AP RL FH, keep empty to use scan size times modelMultiplier
modelParameters.modelMultiplier = 8;
modelParameters.modelFOV = []; % mm, scanner AP RL FH, keep empty to use scan FOV

% Simulate only center of scan region
% Set to region size (in voxels) in MPS (Measurement, Phase, Slice),
% or set to a single scalar n if the region is cubic (equal to [n n n]).
% Set to false to disable.
modelParameters.smallFOVsize = [21 21 21]; % voxels, MPS

% Post processing functions after simulation
options.postProcessKspace = @(x) x .* ringingFilter(size(x)); % Apply ringing filter to simulated k-space
options.postProcessImage = @(x) conj(x); % Convert from left-handed to right-handed spin

%% Load complex image from dicoms
dicomdict('set', 'dicom-dict-philips.txt');

fprintf('Scanning dicom directory...\n');
partitions = readDicomSeries(dicomDir, struct('verbose', true, 'recursive', true));

fprintf('Reading dicom info...\n');
% Read dicominfo
dicomInfo = readDicomSeriesInfo(dicomDir, partitions, struct('SeriesDescription', seriesDescription, 'ImageType', 'ORIGINAL\PRIMARY\R_FFE\R\FFE', 'EchoNumber', echoNumber));

%% Create library
imageParameters = getImageParametersFromDicomPhilips(dicomInfo);
% imageParameters.readoutBandwidth = -imageParameters.readoutBandwidth; % Reverse readout polarity, depends on echo number

library = createLibrary(modelParameters, imageParameters, options);
save -v7.3 'library_Best2301.mat' library

%%
modelParameters.objectFunction = @seedModel_SelectSeed;

library = createLibrary(modelParameters, imageParameters, options);
save -v7.3 'library_I-125_SelectSeed.mat' library
