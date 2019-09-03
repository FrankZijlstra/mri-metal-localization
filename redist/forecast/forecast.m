function [kspaceData] = forecast (model, acquisition, options)
% FORECAST: Fourier-based Off-REsonanCe Artifact simulation in the STeady-state
%
% FORECAST is a fast alternative to Bloch simulation for simulating
% off-resonance effects in steady-state MRI. This function uses a basic
% description of a simulation model and acquisition to perform its
% simulation (see parameter descriptions below). Currently only the model
% and acquistion parameters relevant to off-resonance artifacts in gradient
% echo scans are included, although it is possible to include some other
% parameters (such as T1 and repetition time) yourself (see demo_brainweb.m
% for an example of analytically calculating a steady-state).
%
% See demo_gradientecho.m, demo_spinecho.m, and demo_brainweb.m for some
% usage examples.
%
% If you use FORECAST in your research, please include a reference to our
% MRM paper and a link to the most recent code:
%   F. Zijlstra, J.G. Bouwman, I. Braskute, M.A. Viergever, and P.R.
%   Seevinck, "Fast Fourier-based simulation of off-resonance artifacts in
%   steady-state gradient echo MRI applied to metal object localization",
%   Magn. Reson. Med., 2016
%
% If you have any questions, suggestions, or find any bugs, feel free to
% contact us.
% Frank Zijlstra (f.zijlstra@gmail.com) and Job Bouwman (jgbouwman@hotmail.com), 2017
%
%
% Input parameters:
%   model: Structure with fields defining the simulated model:
%       protonDensity [default: scalar 1]:
%           Defines the effective proton density in arbitrary units. Can be
%           either a scalar constant or a matrix defining a different
%           proton density for each isochromat.
%       deltaB0 [default: scalar 0]:
%           Defines the B0 field shift in Tesla. Can be either a scalar
%           constant or a matrix defining a different field shift for each
%           isochromat. Use calculateFieldShift to calculate this from a
%           susceptibility distribution.
%       T2 [default: scalar Inf]:
%           Defines the T2 values in milliseconds. Can be either a scalar
%           constant or a matrix defining a T2 value for each isochromat.
%           When set to Inf no T2 decay will be simulated.
%       resolution:
%           Vector [X Y Z] containing the resolution (voxel size in mm) for
%           each dimension. Note that Matlab matrices are indexed as
%           M(Y,X,Z).
%   
%   Notes:
%   - At least one of protonDensity, deltaB0 or T2 must be set to a matrix
%   to define the size of the simulation model.
%   - The field of view of the model (resolution * matrix size in XYZ) does
%   not need to match the acquired field of view, but should match exactly
%   with an integer number of voxels in the scan. For example, for a
%   100x100 acquisition matrix with 1 mm resolution, valid model FOVs
%   include 10x10 mm (e.g. 20x20 model with 0.5 mm resolution), 20x80 mm,
%   etc.
%   - The centre of the model matrix is always aligned with the centre of
%   the simulated scan.
%           
%   acquisition: Structure with fields defining the acquisition:
%       kspaceSamplingTimes:
%           Matrix which defines the time (in milliseconds) at which each
%           k-space location is sampled. The size of the matrix defines the
%           size of the k-space itself. Simulation is most efficient for
%           linear Cartesian readouts, use calculateCartesianSamplingTimes
%           to calculate these. These sampling times apply to the T2 decay
%           in the simulation.
%       kspaceSamplingTimesRefocused: [optional]
%           Matrix which defines the refocused sampling time (in
%           milliseconds) at which each k-space location is sampled. These
%           sampling times apply to the delta-B0 induced dephasing in the
%           simulation. If this is not supplied it will be equal to
%           kspaceSamplingTimes. For a basic spin echo simulation:
%             kspaceSamplingTimesRefocused = kspaceSamplingTimes - echoTime
%       resolution:
%           Vector [X Y Z] containing the acquired resolution (voxel size
%           in mm) for each dimension. Note that Matlab matrices are
%           indexed as M(Y,X,Z).
%
%   options: [Optional] Structure with settings
%       verbose [default: false]:
%           When set to true will enable console output of the method and
%           its progress.
%
%   Notes:
%   - An estimate of the acceleration factor relative to Bloch simulation
%   can be obtained with:
%       numel(acquisition.kspaceSamplingTimes) / numel(unique(acquisition.kspaceSamplingTimes))
%
% Output:
%   kspaceData:
%       The simulated k-space data.
%
% Basic usage example:
%   model.protonDensity = phantom(128);
%   model.deltaB0 = randn(128,128)*1e-6; % Random off-resonance
%   model.resolution = [0.5 0.5 1];
%   acquisition.kspaceSamplingTimes = calculateCartesianSamplingTimes('gradientecho', [64 64 1], 'x', 5, 5);
%   acquisition.resolution = [1 1 1];
%   options.verbose = true;
%   kSim = forecast(model, acquisition, options);
%   imshow(abs(ifftc(kSim)),[]);
%



%% Add utils to search path
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

%% Parameter checks
if (nargin < 3)
    options = struct();
end


%% Model parameters

if (~isstruct(model))
    error('Model parameter must be a structure')
end

% If proton-density is not set, assume uniform unit density
if (~isfield(model, 'protonDensity'))
    model.protonDensity = 1;
end

% If delta-B0 is not set, assume no off resonance
if (~isfield(model, 'deltaB0'))
    model.deltaB0 = 0;
end

% If T2 is not set, assume no T2 decay
if (~isfield(model, 'T2'))
    model.T2 = Inf;
end

% Determine model size
model.size = []; % XYZ

if (~isscalar(model.protonDensity))
    model.size([2 1 3]) = size3(model.protonDensity); % XYZ
elseif (~isscalar(model.deltaB0))
    model.size([2 1 3]) = size3(model.deltaB0); % XYZ
elseif (~isscalar(model.T2))
    model.size([2 1 3]) = size3(model.T2); % XYZ
else
    error('forecast: Model is all scalars, at least one of the fields protonDensity, deltaB0, or T2 must be a matrix');
end

% Check that if protonDensity/deltaB0/T2 are given as matrices that the sizes match
if (~isscalar(model.deltaB0) && any(size3(model.deltaB0) ~= model.size([2 1 3])))
    error('forecast: protonDensity and deltaB0 are given as matrix, but their sizes do not match');
end
if (~isscalar(model.T2) && any(size3(model.T2) ~= model.size([2 1 3])))
    error('forecast: protonDensity and/or deltaB0 and T2 are given as matrix, but their sizes do not match');
end

if (~isfield(model, 'resolution'))
    error('forecast: Model structure is missing the resolution field');
end

model.fieldOfView = model.size .* model.resolution; % XYZ


%% Acquisition parameters

if (~isstruct(acquisition))
    error('Acquisition parameter must be a structure')
end

if (~isfield(acquisition, 'kspaceSamplingTimes'))
    error('forecast: Acquisition structure is missing the kspaceSamplingTimes field');
end

if (~isfield(acquisition, 'kspaceSamplingTimesRefocused'))
    acquisition.kspaceSamplingTimesRefocused = acquisition.kspaceSamplingTimes;
end

if (~isequal(size(acquisition.kspaceSamplingTimes), size(acquisition.kspaceSamplingTimesRefocused)))
    error('forecast: The size of kspaceSamplingTimes must be equal to the size kspaceSamplingTimesRefocused');
end

if (~isfield(acquisition, 'resolution'))
    error('forecast: Acquisition structure is missing the resolution field');
end

if (any(model.resolution > acquisition.resolution))
    error('forecast: Model resolution must be higher than acquisition resolution')
end

acquisition.size([2 1 3]) = size3(acquisition.kspaceSamplingTimes); % XYZ
acquisition.fieldOfView = acquisition.size .* acquisition.resolution; % XYZ


%% Options parameters

if (~isstruct(options))
    error('Options parameter must be a structure');
end

if (~isfield(options, 'verbose'))
    options.verbose = false;
end

if (~isfield(options, 'forceGeneral'))
    options.forceGeneral = false;
end

% Set gyromagnetic ratio constant
gm = 42.57747892e6; % Hz / T


%% Check if the model field of view exactly aligns with an integer number of voxels in the scan

scale = model.resolution ./ acquisition.resolution; % XYZ

if (norm(round(model.size .* scale) - model.size .* scale) > 1e-6)
    warning('Model field of view does not exactly match an integer number of voxels in the acquisition matrix. The field of view may be slightly altered.');
end


%%

% Special case: Linear cartesian readout
% Check whether k-space sampling times describe a linear cartesian readout
% in either the X, Y, or Z dimension. A special case simulation will be
% slightly more efficient in these cases.
tmp1 = repmat(acquisition.kspaceSamplingTimes(:,1,1),[1 acquisition.size(1) acquisition.size(3)]);
tmp2 = repmat(acquisition.kspaceSamplingTimesRefocused(:,1,1),[1 acquisition.size(1) acquisition.size(3)]);
useCartesianY = all(acquisition.kspaceSamplingTimes(:) == tmp1(:)) && all(acquisition.kspaceSamplingTimesRefocused(:) == tmp2(:));
tmp1 = repmat(acquisition.kspaceSamplingTimes(1,:,1),[acquisition.size(2) 1 acquisition.size(3)]);
tmp2 = repmat(acquisition.kspaceSamplingTimesRefocused(1,:,1),[acquisition.size(2) 1 acquisition.size(3)]);
useCartesianX = all(acquisition.kspaceSamplingTimes(:) == tmp1(:)) && all(acquisition.kspaceSamplingTimesRefocused(:) == tmp2(:));
tmp1 = repmat(acquisition.kspaceSamplingTimes(1,1,:),[acquisition.size(2) acquisition.size(1) 1]);
tmp2 = repmat(acquisition.kspaceSamplingTimesRefocused(1,1,:),[acquisition.size(2) acquisition.size(1) 1]);
useCartesianZ = all(acquisition.kspaceSamplingTimes(:) == tmp1(:)) && all(acquisition.kspaceSamplingTimesRefocused(:) == tmp2(:));


% Special case: K-Space sampling times are all equal
% Simulation will be more efficient with the generalized simulation, even
% in the linear Cartesian case.
if (all(acquisition.kspaceSamplingTimes(:) == acquisition.kspaceSamplingTimes(1)) && all(acquisition.kspaceSamplingTimesRefocused(:) == acquisition.kspaceSamplingTimesRefocused(1)))
    if (options.verbose)
        fprintf('Using special case: Constant sampling time\n');
    end
    useCartesianY = false;
    useCartesianX = false;
    useCartesianZ = false;
end

% Special case: No field inhomogeneities and no T2
if (all(model.deltaB0(:)==0) && all(model.T2 == Inf))
    if (options.verbose)
        fprintf('Using special case: No field inhomogeneities and no T2\n');
    end
    useCartesianY = false;
    useCartesianX = false;
    useCartesianZ = false;
    acquisition.kspaceSamplingTimes = zeros(size(acquisition.kspaceSamplingTimes));
    acquisition.kspaceSamplingTimesRefocused = zeros(size(acquisition.kspaceSamplingTimesRefocused));
end

% User forced the general simulation
if (options.forceGeneral)
    useCartesianY = false;
    useCartesianX = false;
    useCartesianZ = false;
end


%% Special case: Cartesian Y (note, Y is 1st dimension in Matlab matrices)
if (useCartesianY)
    if (options.verbose)
        fprintf('Using special case: Cartesian Y\n');
    end
    
    nReadout = acquisition.size(2); % Number of readout samples
    
    % Calculate coefficients for DFT
    fc = fftCenter(model.size(2));
    readoutCoeff = (-fc+1:(-fc + model.size(2))).' / model.size(2);
    readoutCoeff = readoutCoeff * model.fieldOfView(2) / acquisition.fieldOfView(2);
    
    % Allocate matrix from intermediate result after DFT
    kSimulation = zeros([nReadout model.size(1) model.size(3)]); % YXZ
    
    % Perform frequency encoding
    for I=1:nReadout
        t_B0 = acquisition.kspaceSamplingTimesRefocused(I,1,1);
        t_T2 = acquisition.kspaceSamplingTimes(I,1,1);
        
        if (options.verbose)
            fprintf('Simulating time %d / %d: %f ms / %f ms\n', I, nReadout, t_B0, t_T2);
        end
        
        % Analytical description of transverse magnetization without encoding at time t
        Mtrans = model.protonDensity .* exp(-2i * pi * gm * model.deltaB0 * (t_B0/1000) - t_T2 ./ model.T2);
        
        % Single frequency DFT in Y
        readEncoding = exp(-2i * pi * readoutCoeff * (I-fftCenter(nReadout)));
        kSimulation(I,:,:) = sum(Mtrans .* repmat(readEncoding, [1 model.size(1) model.size(3)]),1);
    end

    % Perform phase encoding with FFTs
    kSimulation = fft1c(fft1c(kSimulation,[],2),[],3);
    
    % Crop k-space to lower the resolution to the scan resolution
    kspaceData = crop(kSimulation, [nReadout round(model.size(1)*scale(1)) round(model.size(3)*scale(3))]);

    % Correct intensity
    intensityScale = length(readEncoding) / sqrt(nReadout) * sqrt(size(kSimulation,2) / round(model.size(1)*scale(1))) * sqrt(size(kSimulation,3) / round(model.size(3)*scale(3)));
    kspaceData = kspaceData / intensityScale;
    
    % Correct field of view by either cropping or zero-padding in image-space
    kspaceData = fftc(zeroPadOrCrop(ifftc(kspaceData), acquisition.size([2 1 3])));  
    return
end

%% Special case: Cartesian X (note, X is 2nd dimension in Matlab matrices)
if (useCartesianX)
    if (options.verbose)
        fprintf('Using special case: Cartesian X\n');
    end
       
    nReadout = acquisition.size(1); % Number of readout samples

    % Calculate coefficients for DFT
    fc = fftCenter(model.size(1));
    readoutCoeff = (-fc+1:(-fc + model.size(1))) / model.size(1);
    readoutCoeff = readoutCoeff * model.fieldOfView(1) / acquisition.fieldOfView(1);
    
    % Allocate matrix from intermediate result after DFT
    kSimulation = zeros([model.size(2) nReadout model.size(3)]); % YXZ
    
    % Perform frequency encoding
    for I=1:nReadout
        t_B0 = acquisition.kspaceSamplingTimesRefocused(1,I,1);
        t_T2 = acquisition.kspaceSamplingTimes(1,I,1);
        
        if (options.verbose)
            fprintf('Simulating time %d / %d: %f ms / %f ms\n', I, nReadout, t_B0, t_T2);
        end

        % Analytical description of transverse magnetization without encoding at time t
        Mtrans = model.protonDensity .* exp(-2i * pi * gm * model.deltaB0 * (t_B0/1000) - t_T2 ./ model.T2);
        
        % Single frequency DFT in X
        readEncoding = exp(-2i * pi * readoutCoeff * (I-fftCenter(nReadout)));       
        kSimulation(:,I,:) = sum(Mtrans .* repmat(readEncoding, [model.size(2) 1 model.size(3)]),2);
    end

    % Perform phase encoding with FFTs
    kSimulation = fft1c(fft1c(kSimulation,[],1),[],3);
    
    % Crop k-space to lower the resolution to the scan resolution
    kspaceData = crop(kSimulation, [round(model.size(2)*scale(2)) nReadout round(model.size(3)*scale(3))]);

    % Correct intensity
    intensityScale = length(readEncoding) / sqrt(nReadout) * sqrt(size(kSimulation,1) / round(model.size(2)*scale(2))) * sqrt(size(kSimulation,3) / round(model.size(3)*scale(3)));
    kspaceData = kspaceData / intensityScale;
    
    % Correct field of view by either cropping or zero-padding in image-space
    kspaceData = fftc(zeroPadOrCrop(ifftc(kspaceData), acquisition.size([2 1 3])));  
    return
end

%% Special case, cartesian Z
if (useCartesianZ)
    if (options.verbose)
        fprintf('Using special case: Cartesian Z\n');
    end
    
    nReadout = acquisition.size(3); % Number of readout samples
    
    % Calculate coefficients for DFT
    fc = fftCenter(model.size(3));
    readoutCoeff = (-fc+1:(-fc + model.size(3))) / model.size(3);
    readoutCoeff = readoutCoeff * model.fieldOfView(3) / acquisition.fieldOfView(3);
    readoutCoeff = permute(readoutCoeff,[1 3 2]);
    
    % Allocate matrix from intermediate result after DFT
    kSimulation = zeros([model.size(2) model.size(1) nReadout]); % YXZ
    
    % Perform frequency encoding
    for I=1:nReadout
        t_B0 = acquisition.kspaceSamplingTimesRefocused(1,1,I);
        t_T2 = acquisition.kspaceSamplingTimes(1,1,I);
        
        if (options.verbose)
            fprintf('Simulating time %d / %d: %f ms / %f ms\n', I, nReadout, t_B0, t_T2);
        end

        % Analytical description of transverse magnetization without encoding at time t
        Mtrans = model.protonDensity .* exp(-2i * pi * gm * model.deltaB0 * (t_B0/1000) - t_T2 ./ model.T2);
        
        % Single frequency DFT in Z
        readEncoding = exp(-2i * pi * readoutCoeff * (I-fftCenter(nReadout)));
        kSimulation(:,:,I) = sum(Mtrans .* repmat(readEncoding, [model.size(2) model.size(1) 1]),3);
    end

    % Perform phase encoding with FFTs
    kSimulation = fft1c(fft1c(kSimulation,[],1),[],2);
    
    % Crop k-space to lower the resolution to the scan resolution
    kspaceData = crop(kSimulation, [round(model.size(2)*scale(2)) round(model.size(1)*scale(1)) nReadout]);

    % Correct intensity
    intensityScale = length(readEncoding) / sqrt(nReadout) * sqrt(size(kSimulation,1) / round(model.size(2)*scale(2))) * sqrt(size(kSimulation,2) / round(model.size(1)*scale(1)));
    kspaceData = kspaceData / intensityScale;
    
    % Correct field of view by either cropping or zero-padding in image-space
    kspaceData = fftc(zeroPadOrCrop(ifftc(kspaceData), acquisition.size([2 1 3])));  
    return
end


%% Generalized Cartesian simulation
if (options.verbose)
    fprintf('Using generalized Cartesian simulation\n');
end

timeTuples = [acquisition.kspaceSamplingTimesRefocused(:) acquisition.kspaceSamplingTimes(:)];
[Ts,~,Tsind] = unique(timeTuples, 'rows');

kspaceData = zeros(acquisition.size([2 1 3]));
for I=1:size(Ts,1)
    t_B0 = Ts(I,1);
    t_T2 = Ts(I,2);
    
    if (options.verbose)
        fprintf('Simulating time %d / %d: %f ms / %f ms\n', I, size(Ts,1), t_B0, t_T2);
    end
    
    % Analytical description of transverse magnetization without encoding at time t
    Mtrans = model.protonDensity .* exp(-2i * pi * gm * model.deltaB0 * (t_B0/1000) - t_T2 ./ model.T2);
    
    % Spatially encode in all dimensions
    kSimulation = fftc(Mtrans);

    % Crop k-space to match scan resolution
    k = crop(kSimulation, round(model.size([2 1 3]) .* scale([2 1 3])));

    % Correct intensity
    k = k / prod(sqrt(size3(kSimulation) ./ round(model.size([2 1 3]) .* scale([2 1 3]))));

    % Correct field of view by either cropping or zero-padding in image-space
    k = fftc(zeroPadOrCrop(ifftc(k), acquisition.size([2 1 3])));

    % Select elements from k-space that were encoded at time t and save
    % them in the final simulated k-space
    kspaceData(Tsind==I) = k(Tsind==I);
end


end

% Size function that ensures the returned number of dimensions is at least 3
function D = size3 (M)

D = size(M);
if (length(D) < 3)
    D(3) = 1;
end

end
