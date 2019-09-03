function [kspaceSamplingTimes, kspaceSamplingTimesRefocused] = calculateCartesianSamplingTimes (method, matrixSize, readoutDimension, echoTime, readoutInterval)

% calculateCartesianSamplingTimes:
% Function to calculate the time after excitation when k-space locations
% are sampled for linear Cartesian readouts.
%
% Input:
% method: 'gradientecho' or 'ge' for gradient echo sampling times:
%            kspaceSamplingTimesRefocused == kspaceSamplingTimes
%         'spinecho' or 'se' for spin echo sampling times:
%            kspaceSamplingTimesRefocused == kspaceSamplingTimes - echoTime
% matrixSize: [X Y Z] size of the k-space matrix
% readoutDimension: 'x' or 'y' or 'z' dimension along which the readout is performed
% echoTime: Echo time in milliseconds
% readoutInterval: Length of the readout interval in milliseconds
%
% Output:
% kspaceSamplingTimes: Matrix of sampling times per k-space location
%     (applies to T2 decay)
% kspaceSamplingTimesRefocused: Matrix of refocused sampling times per
%     k-space location (applies to delta-B0 induced dephasing)
%
% Frank Zijlstra (f.zijlstra@gmail.com) and Job Bouwman (jgbouwman@hotmail.com), 2017

%% Add utils to search path
addpath(fullfile(fileparts(mfilename('fullpath')), 'utils'));

if (strcmpi(method, 'gradientecho') || strcmpi(method, 'ge'))

    if (strcmpi(readoutDimension, 'x'))
        nReadout = matrixSize(1);
        fc = fftCenter(nReadout);
        kspaceSamplingTimes = repmat((-fc+1:-fc+nReadout), [matrixSize(2) 1 matrixSize(3)]);
    elseif (strcmpi(readoutDimension, 'y'))
        nReadout = matrixSize(2);
        fc = fftCenter(nReadout);
        kspaceSamplingTimes = repmat(permute((-fc+1:-fc+nReadout),[2 1]), [1 matrixSize(1) matrixSize(3)]);
    elseif (strcmpi(readoutDimension, 'z'))
        nReadout = matrixSize(3);
        fc = fftCenter(nReadout);
        kspaceSamplingTimes = repmat(permute((-fc+1:-fc+nReadout),[1 3 2]), [matrixSize(2) matrixSize(1) 1]);
    else
        error('Invalid dimension');
    end

    dT = readoutInterval / nReadout;
    kspaceSamplingTimes = kspaceSamplingTimes * dT + echoTime;
    kspaceSamplingTimesRefocused = kspaceSamplingTimes;
    
elseif (strcmpi(method, 'spinecho') || strcmpi(method, 'se'))
    
    if (strcmpi(readoutDimension, 'x'))
        nReadout = matrixSize(1);
        fc = fftCenter(nReadout);
        kspaceSamplingTimes = repmat((-fc+1:-fc+nReadout), [matrixSize(2) 1 matrixSize(3)]);
    elseif (strcmpi(readoutDimension, 'y'))
        nReadout = matrixSize(2);
        fc = fftCenter(nReadout);
        kspaceSamplingTimes = repmat(permute((-fc+1:-fc+nReadout),[2 1]), [1 matrixSize(1) matrixSize(3)]);
    elseif (strcmpi(readoutDimension, 'z'))
        nReadout = matrixSize(3);
        fc = fftCenter(nReadout);
        kspaceSamplingTimes = repmat(permute((-fc+1:-fc+nReadout),[1 3 2]), [matrixSize(2) matrixSize(1) 1]);
    else
        error('Invalid dimension');
    end

    dT = readoutInterval / nReadout;
    kspaceSamplingTimesRefocused = kspaceSamplingTimes * dT;
    kspaceSamplingTimes = kspaceSamplingTimes * dT + echoTime;
    
else
    error('Unknown acquisition method');
end

end
