clear all
close all

dicomDir = '../testdata/cylinder/';
seriesDescription = 'FFE 1mm 0 deg wfs1';

circularShutter = true;

%% Add necessary paths
addpath('../utils')
addpath('../redist/dicomseries-matlab/dicomseries')
addpath('../simulation')
addpath('../localization')
addpath('../redist/forecast/utils')

%% Load library
load('library.mat');

%% Load complex image from dicoms
dicomdict('set', 'dicom-dict-philips.txt');

fprintf('Scanning dicom directory...\n');
partitions = readDicomSeries(dicomDir, struct('verbose', true, 'recursive', true));

fprintf('Reading images...\n');
% Read echo number 2 with specified protocol name and complex image types
[realImage, realInfo] = readDicomSeriesImage(dicomDir, partitions, struct('SeriesDescription', seriesDescription, 'ImageType', 'ORIGINAL\PRIMARY\R_FFE\R\FFE', 'EchoNumber', 1));
[imagImage, imagInfo] = readDicomSeriesImage(dicomDir, partitions, struct('SeriesDescription', seriesDescription, 'ImageType', 'ORIGINAL\PRIMARY\I_FFE\I\FFE', 'EchoNumber', 1));

fprintf('Rescaling images...\n');
realImage = rescaleDicomImage(realImage, realInfo);
imagImage = rescaleDicomImage(imagImage, imagInfo);

dicomImage = realImage + 1i * imagImage;
clear realImage imagImage

dicomImage = dicomImage / max(abs(dicomImage(:)));

dicomInfo = realInfo; % For accessing general dicom tags
worldMatrix = getWorldMatrixFromDicom(dicomInfo);



%% Template matching
data = struct();
data.library = library;
data.image = dicomImage;

if (circularShutter)
    kSpaceMask = createCircularShutter(size(data.image), dicomInfo);
else
    kSpaceMask = ones(size(data.image));
end

fprintf('Matching library...\n');
[data.match, data.templateId] = matchLibrary(data.image, data.library, kSpaceMask);

[~,ind] = max(data.match(:));
[y,x,z] = ind2sub(size(data.match), ind);
position = worldMatrix * [x-1, y-1, z-1, 1].';
orientation = library.templates(data.templateId(ind)).orientationVector;

%% Visualization
close all
SLICE = 29;

figure
hold on
orthoplotWorldMatrix(abs(data.image) / max(abs(data.image(:))), fftCenter(size(data.image,1)), fftCenter(size(data.image,2)), SLICE, worldMatrix)

msh = cylinderMeshWorld(25/2, 80, position(1:3).', orientation);
p = patch(msh);
set(p,'FaceColor','red','EdgeColor','none');
axis equal
campos([-2000 1800 1500])
camzoom(2)
