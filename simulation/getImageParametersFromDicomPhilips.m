function [imageParameters] = getImageParametersFromDicomPhilips (dicomInfo)

imageParameters = struct();
IOP = getDicomAttribute(dicomInfo, 'ImageOrientationPatient');
patientPosition = getDicomAttribute(dicomInfo, 'PatientPosition');

phaseEncodingDirection = getDicomAttribute(dicomInfo, 'InPlanePhaseEncodingDirection');
voxelSize = [getDicomAttribute(dicomInfo, 'PixelSpacing'); getDicomAttribute(dicomInfo, 'SliceThickness')];

if strcmp(phaseEncodingDirection, 'ROW')
    imageParameters.readoutDimension = 'y';
else
    imageParameters.readoutDimension = 'x';
end

% patient == scanner only in case of patientPosition == HFS
imageParameters.readoutOrientation = IOP([2 1 3]); % AP RL FH (patient)
imageParameters.phaseOrientation = IOP([5 4 6]); % AP RL FH (patient)

imageParameters.imageSize = double([getDicomAttribute(dicomInfo, 'Columns') getDicomAttribute(dicomInfo, 'Rows') getDicomAttribute(dicomInfo, 'MRSeriesNrOfSlices')]); % Y X Z
imageParameters.imageFOV = imageParameters.imageSize .* voxelSize.'; % Y X Z

imageParameters.readoutBandwidth = getDicomAttribute(dicomInfo, 'PixelBandwidth'); % Hz/Pixel (add minus for opposite readout polarity)
imageParameters.echoTime = getDicomAttribute(dicomInfo, 'EchoTime');
imageParameters.scanType = 'gradientecho'; % or 'spinecho'
imageParameters.fieldStrength = getDicomAttribute(dicomInfo, 'MRSeriesMagneticFieldStrength'); % T

end
