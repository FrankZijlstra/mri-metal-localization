function [M,R] = getWorldMatrixFromDicom (info)
%This function calculates the 4x4 transform matrix from the image
%coordinates to patient coordinates. 

ipp = getDicomAttribute(info, 'ImagePositionPatient');
iop = getDicomAttribute(info, 'ImageOrientationPatient');
ps = getDicomAttribute(info, 'PixelSpacing');
ss = getDicomAttribute(info, 'SpacingBetweenSlices');
if isempty(ss)
    ss = getDicomAttribute(info, 'SliceThickness');
end

Tipp = [1 0 0 ipp(1);
        0 1 0 ipp(2);
        0 0 1 ipp(3);
        0 0 0 1];

r = iop(1:3);
c = iop(4:6);
s = cross(r', c');

R = [r(1) c(1) s(1) 0;
     r(2) c(2) s(2) 0;
     r(3) c(3) s(3) 0;
     0 0 0 1];

S = [ps(2) 0 0 0;
     0 ps(1) 0 0;
     0 0 abs(ss) 0;
     0 0 0 1];

T0 = [1 0 0 0;
      0 1 0 0;
      0 0 1 0;
      0 0 0 1];
  
M = Tipp * R * S * T0;

end