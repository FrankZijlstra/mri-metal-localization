clear all
close all

dicomDir = '../testdata/brachyseeds/';
seriesDescription = 'mFFE';
echoNumber = 2;

numSeeds = 3;
numberOfCandidates = 1000;
circularShutter = true;

% Bounding box defining a region of interest. This can also be based on a
% segmentation, etc. Note that you do need to keep a boundary of half the
% template size, since no template matching will be performed in the
% boundary.
% The bounding box can be set to [1 1 1], size(dicomImage) for full image
% (or by removing all bounding box code below)
boundingBoxMin = [1 1 1]; % YXZ
boundingBoxMax = [192 192 64]; % YXZ


%% Add necessary paths
addpath('../utils')
addpath('../redist/dicomseries-matlab/dicomseries')
addpath('../simulation')
addpath('../localization')
addpath('../redist/forecast/utils')

%% Load library
% load('library_Best2301.mat');
load('library_I-125_SelectSeed.mat');

%% Load complex image from dicoms
dicomdict('set', 'dicom-dict-philips.txt');

fprintf('Scanning dicom directory...\n');

partitions = readDicomSeries(dicomDir, struct('verbose', true, 'recursive', true));

fprintf('Reading images...\n');
% Read echo number 2 with specified protocol name and complex image types
[realImage, dicomInfo] = readDicomSeriesImage(dicomDir, partitions, struct('SeriesDescription', seriesDescription, 'ImageType', 'ORIGINAL\PRIMARY\R_FFE\R\FFE', 'EchoNumber', echoNumber));
[imagImage, imagInfo] = readDicomSeriesImage(dicomDir, partitions, struct('SeriesDescription', seriesDescription, 'ImageType', 'ORIGINAL\PRIMARY\I_FFE\I\FFE', 'EchoNumber', echoNumber));

fprintf('Rescaling images...\n');
realImage = rescaleDicomImage(realImage, dicomInfo);
imagImage = rescaleDicomImage(imagImage, imagInfo);

dicomImage = realImage + 1i * imagImage;
clear realImage imagImage

dicomImage = dicomImage / max(abs(dicomImage(:)));
worldMatrix = getWorldMatrixFromDicom(dicomInfo);


%% Template matching
data = struct();
data.library = library;

% Crop image to bounding box
data.image = dicomImage(boundingBoxMin(1):boundingBoxMax(1),boundingBoxMin(2):boundingBoxMax(2),boundingBoxMin(3):boundingBoxMax(3));

% Calculate world matrix for cropped image
shiftMatrix = eye(4);
shiftMatrix(1:3,4) = boundingBoxMin-1;
data.worldMatrix = worldMatrix * shiftMatrix;

if (circularShutter)
    kSpaceMask = createCircularShutter(size(data.image), dicomInfo);
else
    kSpaceMask = ones(size(data.image));
end

fprintf('Matching library...\n');
[data.match, data.templateId] = matchLibrary(data.image, data.library, kSpaceMask);
data.match = removeInvalidMatches(data.match, data.library);


%% Linear regression on candidate dectections
tmp = sort(abs(data.match(:)),'descend');
phaseCorrelationThreshold = tmp(numberOfCandidates);

data.mask = abs(data.match) >= phaseCorrelationThreshold;

fprintf('Calculating regression errors...\n');
[data.regressionError, data.regressionSlope, data.regressionOffset] = matchLibraryRegression(data.image, data.templateId, data.library, data.mask);


%% Detect seeds (note: can change data)
seedPositions = [];
seedOrientations = [];
K = 0;

while (K < numSeeds)
    [data, position, orientation] = detectSeed(data);
    if (isempty(position))
        break;
    end
    K = K + 1;
    fprintf('Detected seed %d\n', K);
    seedPositions(K,:) = position; %#ok<SAGROW>
    seedOrientations(K,:) = orientation; %#ok<SAGROW>
end


%% Visualization
close all
SLICE = 23;

figure
hold on
orthoplotWorldMatrix(abs(data.image) / max(abs(data.image(:))), fftCenter(size(data.image,1)), fftCenter(size(data.image,2)), SLICE, data.worldMatrix)

for I=1:size(seedPositions,1)
    msh = cylinderMeshWorld(0.4, 4, seedPositions(I,:), seedOrientations(I,:) .* [1 1 -1]);
    p = patch(msh);
    set(p,'FaceColor','red','EdgeColor','none');
end
axis equal
campos([-200 300 140])
camzoom(4)