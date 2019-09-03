function [library] = createLibrary (modelParameters, imageParameters, options)



% TODO: The sign of readout/phase orientation is not used
% TODO: Angulated orientations not supported

% TODO: Allow modelMultiplier to be a 3 element vector?
% TODO: Allow specification of model voxel size? (i.e. you provide the
% voxel size and it determines the model size automatically)

% TODO: Allow specification of FOV/size in MPS, instead of AP/RL/FH?
% (options.coordinateSystem = 'MPS' or 'scanner'?)


smallFOVsize = modelParameters.smallFOVsize; % voxels, MPS
simulateSmallFOV = ~isempty(smallFOVsize);


% Post processing functions after simulation
if (~isfield(options, 'postProcessKspace'))
    options.postProcessKspace = @(x) x;
end
if (~isfield(options, 'postProcessImage'))
    options.postProcessImage = @(x) x;
end

postProcessKspace = options.postProcessKspace;
postProcessImage = options.postProcessImage;

%% Initialize library with parameters
library = struct();

library.info.modelSize = modelParameters.modelSize;
library.info.modelFOV = modelParameters.modelFOV;
library.info.modelMultiplier = modelParameters.modelMultiplier;
library.info.simulateSmallFOV = simulateSmallFOV;
library.info.smallFOVsize = modelParameters.smallFOVsize;

library.info.postProcessKspace = postProcessKspace;
library.info.postProcessImage = postProcessImage;

library.info.B0 = imageParameters.fieldStrength;


%%

FOVPE = imageParameters.imageFOV(1); % mm
FOVreadout = imageParameters.imageFOV(2); % mm
FOVslice = imageParameters.imageFOV(3); % mm

NPE = imageParameters.imageSize(1);
Nreadout = imageParameters.imageSize(2);
Nslices = imageParameters.imageSize(3);

scanVoxelSize = [FOVreadout/Nreadout FOVPE/NPE FOVslice/Nslices]; % MPS

library.info.scanFOV = [FOVreadout FOVPE FOVslice];
library.info.scanSize = [Nreadout NPE Nslices];
library.info.scanResolution = scanVoxelSize;

% Change scan size and FOV if we simulate a smaller FOV
if (simulateSmallFOV)
    if (isscalar(smallFOVsize))
        NPE = smallFOVsize;
        Nreadout = smallFOVsize;
        Nslices = smallFOVsize;
    else
        NPE = smallFOVsize(2);
        Nreadout = smallFOVsize(1);
        Nslices = smallFOVsize(3);
    end
    
    FOVreadout = Nreadout * scanVoxelSize(1);
    FOVPE = NPE * scanVoxelSize(2);
    FOVslice = Nslices * scanVoxelSize(3);
    
    scanVoxelSize = [FOVreadout/Nreadout FOVPE/NPE FOVslice/Nslices]; % MPS
end

readoutDir = find(imageParameters.readoutOrientation);
phaseDir = find(imageParameters.phaseOrientation);
sliceDir = find(cross(imageParameters.phaseOrientation, imageParameters.readoutOrientation));

scanPermutation = [phaseDir readoutDir sliceDir]; % Value tells which direction each image dimension in PMS is (1 = AP, 2 = RL, 3 = FH)
scanPermutationInv = [find(scanPermutation == 1) find(scanPermutation == 2) find(scanPermutation == 3)]; % From AP RL FH to phase readout slice


library.info.scanPermutation = scanPermutation; % PMS (1 = AP, 2 = RL, 3 = FH)
library.info.scanPermutationInv = scanPermutationInv; % AP RL FH (1 = P, 2 = M, 3 = S)


%% Calculate sampling times

echoTime = imageParameters.echoTime; % ms
dT = 1000/imageParameters.readoutBandwidth;

acquisition = struct();
[acquisition.kspaceSamplingTimes, acquisition.kspaceSamplingTimesRefocused] = calculateCartesianSamplingTimes(imageParameters.scanType, [Nreadout NPE Nslices], imageParameters.readoutDimension, echoTime, dT);
acquisition.resolution = library.info.scanResolution;

library.info.echoTime = echoTime;


%%
modelSize = modelParameters.modelSize;
modelMultiplier = modelParameters.modelMultiplier;

if (isempty(modelSize))
    modelSize(readoutDir) = Nreadout * modelMultiplier;
    modelSize(phaseDir) = NPE * modelMultiplier;
    modelSize(sliceDir) = Nslices * modelMultiplier;
    
    fprintf('Automatically set model size: AP: %d RL: %d FH: %d\n', modelSize(1), modelSize(2), modelSize(3));
end

Ny = modelSize(1); % AP
Nx = modelSize(2); % RL
Nz = modelSize(3); % FH

% Object will be AP RL FH(B0)

scanFOVAPRLFH(readoutDir) = FOVreadout;
scanFOVAPRLFH(phaseDir) = FOVPE;
scanFOVAPRLFH(sliceDir) = FOVslice;

modelFOV = modelParameters.modelFOV;

if (~isempty(modelFOV))
    FOVy = modelFOV(1); % AP
    FOVx = modelFOV(2); % RL
    FOVz = modelFOV(3); % FH
else
    FOVy = scanFOVAPRLFH(1); % AP
    FOVx = scanFOVAPRLFH(2); % RL
    FOVz = scanFOVAPRLFH(3); % FH
    
    fprintf('Automatically set model FOV: AP: %.2f RL: %.2f FH: %.2f\n', FOVy, FOVx, FOVz);
end

FOVy = [-FOVy/2 FOVy/2]; % AP
FOVx = [-FOVx/2 FOVx/2]; % RL
FOVz = [-FOVz/2 FOVz/2]; % FH

XYZtoYXZ = [2 1 3];

scanVoxelSizeAPRLFH = scanVoxelSize(XYZtoYXZ(scanPermutationInv));

% Make FOV a multiple of scanVoxelSize
% TODO: Give warning if this changes the FOV?

FOVAPRLFH = [FOVy(2)-FOVy(1) FOVx(2)-FOVx(1) FOVz(2)-FOVz(1)]; % AP RL FH
simResolution = FOVAPRLFH ./ modelSize; % AP RL FH
scale = simResolution ./ scanVoxelSizeAPRLFH; % AP RL FH

FOV = round(scale.*modelSize).*scanVoxelSizeAPRLFH;

FOVx = [-FOV(2)/2 FOV(2)/2];
FOVy = [-FOV(1)/2 FOV(1)/2];
FOVz = [-FOV(3)/2 FOV(3)/2];

FOVAPRLFH = [FOVy(2)-FOVy(1) FOVx(2)-FOVx(1) FOVz(2)-FOVz(1)]; % AP RL FH

%%

voxelSize = [(FOVy(2)-FOVy(1))/Ny (FOVx(2)-FOVx(1))/Nx (FOVz(2)-FOVz(1))/Nz]; % AP RL FH

xSpace = calculateSimulationCoordinates(FOVx, Nx);
ySpace = calculateSimulationCoordinates(FOVy, Ny);
zSpace = calculateSimulationCoordinates(FOVz, Nz);

fprintf('Effective simulation FOV:\nx: %f %f\ny: %f %f\nz: %f %f\n', xSpace(1)-voxelSize(2)/2,xSpace(end)+voxelSize(2)/2, ySpace(1)-voxelSize(1)/2,ySpace(end)+voxelSize(1)/2, zSpace(1)-voxelSize(3)/2,zSpace(end)+voxelSize(2)/2);

%%

library.info.templateFOV = [FOVreadout FOVPE FOVslice]; % MPS
library.info.templateSize = [Nreadout NPE Nslices]; % MPS
library.info.templateResolution = scanVoxelSize; % MPS

library.info.simFOV = FOVAPRLFH(scanPermutation(XYZtoYXZ)); % MPS
library.info.simSize = modelSize(scanPermutation(XYZtoYXZ)); % MPS
library.info.simResolution = library.info.simFOV ./ library.info.simSize; % MPS

%%

fprintf('Simulating background... ');
tic
% Generate object model
model = struct();
model.resolution = library.info.simResolution;

[model.protonDensity, model.deltaB0, model.T2] = modelParameters.objectFunction(xSpace, ySpace, zSpace, [], voxelSize, imageParameters.fieldStrength);

% Transform model matrices from AP RL FH to PMS
if (~isscalar(model.protonDensity))
    model.protonDensity = permute(model.protonDensity, scanPermutation);
end
if (~isscalar(model.deltaB0))
    model.deltaB0 = permute(model.deltaB0, scanPermutation);
end
if (~isscalar(model.T2))
    model.T2 = permute(model.T2, scanPermutation);
end

% Run scan simulation
kSim = forecast(model, acquisition, struct('verbose', true));
toc

% Post-processing
if (~isempty(postProcessKspace))
    kSim = postProcessKspace(kSim);
end

bg = ifftc(kSim);

if (~isempty(postProcessImage))
    bg = postProcessImage(bg);
end

library.background = bg;


%%

library.info.numberOfTemplates = size(modelParameters.objectOrientations,1);

for I=1:size(modelParameters.objectOrientations,1)
    fprintf('Template: %d\n', I);
    orientationVector = modelParameters.objectOrientations(I,:);
    
    fprintf('Simulating... ');
    tic
    
    % Calculate rotation relative to orientation in z-direction
    R = vrrotvec2mat(vrrotvec(orientationVector, [0 0 1])); % Rotation from rotated space to scanner space
   
    % Generate object model
    model = struct();
    model.resolution = library.info.simResolution;
    
    [model.protonDensity, model.deltaB0, model.T2] = modelParameters.objectFunction(xSpace, ySpace, zSpace, R, voxelSize, imageParameters.fieldStrength);

    % Transform model matrices from AP RL FH to PMS
    if (~isscalar(model.protonDensity))
        model.protonDensity = permute(model.protonDensity, scanPermutation);
    end
    if (~isscalar(model.deltaB0))
        model.deltaB0 = permute(model.deltaB0, scanPermutation);
    end
    if (~isscalar(model.T2))
        model.T2 = permute(model.T2, scanPermutation);
    end

    % Run scan simulation   
    kSim = forecast(model, acquisition, struct('verbose', true));
    toc
    
    % Post-processing
    if (~isempty(postProcessKspace))
        kSim = postProcessKspace(kSim);
    end
    
    template = ifftc(kSim);
    
    if (~isempty(postProcessImage))
        template = postProcessImage(template);
    end

    library.templates(I).image = template;
    % TODO: Store orientationVector in dicom world space
    library.templates(I).orientationVector = orientationVector; % TODO: orientation + roll is necessary for non cylindrical objects
end

end
